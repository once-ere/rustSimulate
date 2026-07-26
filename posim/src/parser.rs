//! Parser ("yacc/bison" analog): a recursive-descent grammar compiler
//! that turns the token stream into a postfix instruction program for
//! the stack machine in [`crate::vm`].
//!
//! Grammar (EBNF):
//!
//! ```text
//! command  := "NEW" shape [ "AS" (IDENT | STRING) ]
//!                          [ "{" init { "," init } "}" ]
//!                                               (* AS registers a user
//!                                                  name for the object *)
//!           | "SET" path "=" expr
//!           | "GET" path
//!           | "DEL" NUMBER
//!           | "LIST"
//!           | "STEP" expr                       (* advance by dt      *)
//!           | "RUN" expr [ "STEPS" NUMBER ]     (* advance by t, n outs *)
//!           | "METHOD" ( "ADAMS" | "BDF" | "SPRK" IDENT [ NUMBER ] )
//!           | "ENERGY" | "COM" | "MOMENTUM" | "ANGMOM"
//!           | "LAPLACE" NUMBER
//!           | "RESET" | "HELP"
//!           | "SCENE" scenecmd                  (* graphical scene     *)
//!           | "COLLIDE" [ "ON" | "OFF" ]        (* bare: report status *)
//!           | "CONTACTS"                        (* list last contacts  *)
//!           | "LET" IDENT "=" expr              (* session variable    *)
//!           | "FUNCS"                           (* list user functions *)
//!           | "SHOW" IDENT                      (* print a function    *)
//!           | "BOX" [ "OFF" | expr ]            (* rigid bounding box:
//!                                                  expr = inner side
//!                                                  length; bare = status;
//!                                                  OFF removes it      *)
//!           | expr ;                            (* bare expression     *)
//! scenecmd := "CREATE" [ NUMBER ]               (* open window [port]  *)
//!           | "CLOSE"                           (* aka DESTROY         *)
//!           | "TRANSLATE" term term [ term ]    (* camera dx dy [dz]   *)
//!           | "ROTATE" term term                (* camera dyaw dpitch  *)
//!           | "ZOOM" ( "IN" | "OUT" | term )    (* factor > 1 zooms in *)
//!           | "HIDE" [ NUMBER | "ALL" ]         (* default: ALL        *)
//!           | "SHOW" [ NUMBER | "ALL" ]
//!           | "REFRESH"                         (* re-sync from state  *)
//!           | "REDRAW"                          (* re-send full scene  *)
//!           | "START" | "STOP" | "PAUSE" | "REVERSE"
//!           | "RESET"                           (* re-initialize: all
//!                                                  values and the time
//!                                                  return to initial;
//!                                                  START re-starts    *)
//!           | "SET_TIME_STEP" term              (* args are term-level:
//!                                                  -5 is negative five;
//!                                                  parenthesize sums  *)
//!           | "STATUS" | "EVENTS" ;
//! shape    := "POINT" | "SPHERE" | "CUBOID" | "TORUS" | "DISK" | "CYLINDER"
//!           | "DUMBBELL" ;                      (* two solid spheres +
//!                                                  a rigid rod, one
//!                                                  rigid body          *)
//!
//! (* User-defined functions are a LINE FORM handled before this
//!    grammar: DEF name(param [= default], ...) { body } — the body is
//!    newline/;-separated commands using the parameters as variables;
//!    each body line must itself satisfy this grammar. Invocation uses
//!    the ordinary call syntax name(arg, ...); trailing parameters
//!    take their defaults. *)
//! init     := IDENT "=" expr ;
//! path     := IDENT { "." IDENT } ;             (* objN.field[.x|y|z|w],
//!                                                  system.field,
//!                                                  contactK.field,
//!                                                  name.field for
//!                                                  AS-registered names *)
//! expr     := term { ("+" | "-") term } ;
//! term     := unary { ("*" | "/") unary } ;
//! unary    := "-" unary | atom ;
//! atom     := NUMBER | STRING | "[" expr { "," expr } "]" | "(" expr ")"
//!           | IDENT "(" [ expr { "," expr } ] ")"   (* builtin or user
//!                                                      function call   *)
//!           | path
//!           | IDENT ;                           (* parameter / LET var *)
//! ```

use crate::lexer::{tokenize, Keyword, TokKind, Token};
use crate::vm::{Instr, MethodSpec, NameArg, Path, PathRoot, ShapeKind, Value};

pub struct Parser {
    toks: Vec<Token>,
    pos: usize,
}

/// Compiles exactly ONE expression (never a command form) — the entry
/// point for DEF parameter defaults, so a default like `reset` or
/// `del 0` is a definition-time error instead of a command that runs.
pub fn compile_expression(src: &str) -> Result<Vec<Instr>, String> {
    let toks = tokenize(src)?;
    if toks.is_empty() {
        return Err("expected an expression".to_string());
    }
    let mut p = Parser { toks, pos: 0 };
    let mut prog = Vec::new();
    p.expr(&mut prog)?;
    if let Some(t) = p.peek() {
        return Err(format!(
            "parse error at column {}: unexpected {} after the expression",
            t.col, t.kind
        ));
    }
    Ok(prog)
}

/// Compiles one command line into a stack-machine program.
pub fn compile_line(line: &str) -> Result<Vec<Instr>, String> {
    let toks = tokenize(line)?;
    if toks.is_empty() {
        return Ok(Vec::new());
    }
    let mut p = Parser { toks, pos: 0 };
    let prog = p.command()?;
    if let Some(t) = p.peek() {
        return Err(format!("parse error at column {}: unexpected {}", t.col, t.kind));
    }
    Ok(prog)
}

impl Parser {
    fn peek(&self) -> Option<&Token> {
        self.toks.get(self.pos)
    }

    fn next(&mut self) -> Option<Token> {
        let t = self.toks.get(self.pos).cloned();
        if t.is_some() {
            self.pos += 1;
        }
        t
    }

    fn eat(&mut self, want: &TokKind) -> Result<(), String> {
        match self.next() {
            Some(t) if t.kind == *want => Ok(()),
            Some(t) => Err(format!(
                "parse error at column {}: expected {}, found {}",
                t.col, want, t.kind
            )),
            None => Err(format!("parse error: expected {} at end of line", want)),
        }
    }

    fn eat_keyword(&mut self, kw: Keyword) -> bool {
        if matches!(self.peek(), Some(t) if t.kind == TokKind::Keyword(kw)) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn expect_number(&mut self, what: &str) -> Result<f64, String> {
        match self.next() {
            Some(Token { kind: TokKind::Number(n), .. }) => Ok(n),
            Some(t) => Err(format!(
                "parse error at column {}: expected {what}, found {}",
                t.col, t.kind
            )),
            None => Err(format!("parse error: expected {what} at end of line")),
        }
    }

    fn expect_ident(&mut self, what: &str) -> Result<String, String> {
        match self.next() {
            Some(Token { kind: TokKind::Ident(s), .. }) => Ok(s),
            Some(t) => Err(format!(
                "parse error at column {}: expected {what}, found {}",
                t.col, t.kind
            )),
            None => Err(format!("parse error: expected {what} at end of line")),
        }
    }

    /// A field name may collide with a keyword (`momentum`, `energy`,
    /// `method`, ...): after a `.` or inside `NEW { ... }` both are
    /// accepted.
    fn expect_field(&mut self) -> Result<String, String> {
        match self.next() {
            Some(Token { kind: TokKind::Ident(s), .. }) => Ok(s.to_ascii_lowercase()),
            Some(Token { kind: TokKind::Keyword(k), .. }) => Ok(format!("{k:?}").to_ascii_lowercase()),
            Some(t) => Err(format!(
                "parse error at column {}: expected a field name, found {}",
                t.col, t.kind
            )),
            None => Err("parse error: expected a field name at end of line".to_string()),
        }
    }

    fn command(&mut self) -> Result<Vec<Instr>, String> {
        let mut prog = Vec::new();
        let t = self.peek().cloned().expect("nonempty");
        match t.kind {
            TokKind::Keyword(Keyword::New) => {
                self.pos += 1;
                let shape = match self.next() {
                    Some(Token { kind: TokKind::Keyword(Keyword::Point), .. }) => ShapeKind::Point,
                    Some(Token { kind: TokKind::Keyword(Keyword::Sphere), .. }) => ShapeKind::Sphere,
                    Some(Token { kind: TokKind::Keyword(Keyword::Cuboid), .. }) => ShapeKind::Cuboid,
                    Some(Token { kind: TokKind::Keyword(Keyword::Torus), .. }) => ShapeKind::Torus,
                    Some(Token { kind: TokKind::Keyword(Keyword::Disk), .. }) => ShapeKind::Disk,
                    Some(Token { kind: TokKind::Keyword(Keyword::Cylinder), .. }) => {
                        ShapeKind::Cylinder
                    }
                    Some(Token { kind: TokKind::Keyword(Keyword::Dumbbell), .. }) => {
                        ShapeKind::Dumbbell
                    }
                    Some(t) => {
                        return Err(format!(
                            "parse error at column {}: expected POINT, SPHERE, CUBOID, TORUS, \
                             DISK, CYLINDER or DUMBBELL, found {}",
                            t.col, t.kind
                        ))
                    }
                    None => {
                        return Err("parse error: NEW needs a shape (POINT, SPHERE, CUBOID, \
                                    TORUS, DISK, CYLINDER, DUMBBELL)"
                            .into())
                    }
                };
                prog.push(Instr::NewObject(shape));
                /* optional `AS <name>`: a bare identifier (resolved
                 * against function parameters / LET variables holding a
                 * string, else taken literally) or a string literal */
                let name = if self.eat_keyword(Keyword::As) {
                    match self.next() {
                        Some(Token { kind: TokKind::Ident(n), .. }) => {
                            Some(NameArg::Ident(n.to_ascii_lowercase()))
                        }
                        Some(Token { kind: TokKind::Str(st), .. }) => Some(NameArg::Str(st)),
                        Some(t) => {
                            return Err(format!(
                                "parse error at column {}: AS needs a name (identifier or \
                                 string), found {}",
                                t.col, t.kind
                            ))
                        }
                        None => return Err("parse error: AS needs a name".into()),
                    }
                } else {
                    None
                };
                let mut explicit_inertia = false;
                if matches!(self.peek(), Some(t) if t.kind == TokKind::LBrace) {
                    self.pos += 1;
                    loop {
                        let field = self.expect_field()?;
                        if field == "inertia_tensor" || field == "inverse_inertia_tensor" {
                            explicit_inertia = true;
                        }
                        self.eat(&TokKind::Equals)?;
                        self.expr(&mut prog)?;
                        prog.push(Instr::InitField(field));
                        if matches!(self.peek(), Some(t) if t.kind == TokKind::Comma) {
                            self.pos += 1;
                            continue;
                        }
                        break;
                    }
                    self.eat(&TokKind::RBrace)?;
                }
                prog.push(Instr::FinishNew { recompute_inertia: !explicit_inertia, name });
            }
            TokKind::Keyword(Keyword::Set) => {
                self.pos += 1;
                let path = self.path()?;
                self.eat(&TokKind::Equals)?;
                self.expr(&mut prog)?;
                prog.push(Instr::Store(path));
            }
            TokKind::Keyword(Keyword::Get) => {
                self.pos += 1;
                let path = self.path()?;
                prog.push(Instr::Load(path));
            }
            TokKind::Keyword(Keyword::Del) => {
                self.pos += 1;
                let n = self.expect_number("an object index")?;
                prog.push(Instr::Delete(n as usize));
            }
            TokKind::Keyword(Keyword::List) => {
                self.pos += 1;
                prog.push(Instr::ListObjects);
            }
            TokKind::Keyword(Keyword::Step) => {
                self.pos += 1;
                self.expr(&mut prog)?;
                prog.push(Instr::Step);
            }
            TokKind::Keyword(Keyword::Run) => {
                self.pos += 1;
                self.expr(&mut prog)?;
                let steps = if self.eat_keyword(Keyword::Steps) {
                    self.expect_number("a step count")? as usize
                } else {
                    10
                };
                prog.push(Instr::Run { outputs: steps.max(1) });
            }
            TokKind::Keyword(Keyword::Method) => {
                self.pos += 1;
                let spec = match self.next() {
                    Some(Token { kind: TokKind::Keyword(Keyword::Adams), .. }) => MethodSpec::Adams,
                    Some(Token { kind: TokKind::Keyword(Keyword::Bdf), .. }) => MethodSpec::Bdf,
                    Some(Token { kind: TokKind::Keyword(Keyword::Sprk), .. }) => {
                        let raw = self.expect_ident("an SPRK table name")?;
                        let upper = raw.to_ascii_uppercase();
                        let table = if upper.starts_with("ARKODE_") {
                            upper
                        } else {
                            format!("ARKODE_SPRK_{upper}")
                        };
                        let dt = match self.peek() {
                            Some(Token { kind: TokKind::Number(_), .. }) => {
                                self.expect_number("a fixed step dt")?
                            }
                            _ => 0.01,
                        };
                        MethodSpec::Sprk { table, dt }
                    }
                    Some(t) => {
                        return Err(format!(
                            "parse error at column {}: expected ADAMS, BDF or SPRK, found {}",
                            t.col, t.kind
                        ))
                    }
                    None => return Err("parse error: METHOD needs ADAMS, BDF or SPRK".into()),
                };
                prog.push(Instr::SetMethod(spec));
            }
            TokKind::Keyword(Keyword::Energy) => {
                self.pos += 1;
                prog.push(Instr::Energy);
            }
            TokKind::Keyword(Keyword::Com) => {
                self.pos += 1;
                prog.push(Instr::CenterOfMass);
            }
            TokKind::Keyword(Keyword::Momentum) => {
                self.pos += 1;
                prog.push(Instr::TotalMomentum);
            }
            TokKind::Keyword(Keyword::Angmom) => {
                self.pos += 1;
                prog.push(Instr::TotalAngularMomentum);
            }
            TokKind::Keyword(Keyword::Laplace) => {
                self.pos += 1;
                let n = self.expect_number("an object index")?;
                prog.push(Instr::Laplace(n as usize));
            }
            TokKind::Keyword(Keyword::Reset) => {
                self.pos += 1;
                prog.push(Instr::Reset);
            }
            TokKind::Keyword(Keyword::Help) => {
                self.pos += 1;
                prog.push(Instr::Help);
            }
            TokKind::Keyword(Keyword::Scene) => {
                self.pos += 1;
                self.scene_command(&mut prog)?;
            }
            TokKind::Keyword(Keyword::Collide) => {
                self.pos += 1;
                let mode = match self.peek() {
                    Some(Token { kind: TokKind::Keyword(Keyword::On), .. }) => {
                        self.pos += 1;
                        Some(true)
                    }
                    Some(Token { kind: TokKind::Keyword(Keyword::Off), .. }) => {
                        self.pos += 1;
                        Some(false)
                    }
                    _ => None,
                };
                prog.push(Instr::Collide(mode));
            }
            TokKind::Keyword(Keyword::Contacts) => {
                self.pos += 1;
                prog.push(Instr::Contacts);
            }
            TokKind::Keyword(Keyword::Let) => {
                self.pos += 1;
                let name = self.expect_ident("a variable name")?.to_ascii_lowercase();
                if name == "pi" || name == "tau" {
                    return Err(format!(
                        "`{name}` is a built-in constant and cannot be a LET variable"
                    ));
                }
                self.eat(&TokKind::Equals)?;
                self.expr(&mut prog)?;
                prog.push(Instr::StoreGlobal(name));
            }
            TokKind::Keyword(Keyword::Funcs) => {
                self.pos += 1;
                prog.push(Instr::ListFns);
            }
            TokKind::Keyword(Keyword::Show) => {
                self.pos += 1;
                let name = self.expect_ident("a function name")?.to_ascii_lowercase();
                prog.push(Instr::ShowFn(name));
            }
            TokKind::Keyword(Keyword::Box) => {
                self.pos += 1;
                use crate::vm::BoxMode;
                match self.peek() {
                    Some(Token { kind: TokKind::Keyword(Keyword::Off), .. }) => {
                        self.pos += 1;
                        prog.push(Instr::Box(BoxMode::Off));
                    }
                    None => prog.push(Instr::Box(BoxMode::Status)),
                    _ => {
                        self.expr(&mut prog)?;
                        prog.push(Instr::Box(BoxMode::Create));
                    }
                }
            }
            _ => {
                self.expr(&mut prog)?;
            }
        }
        Ok(prog)
    }

    /// `scenecmd` — the sub-command after the `SCENE` keyword.
    fn scene_command(&mut self, prog: &mut Vec<Instr>) -> Result<(), String> {
        use crate::vm::SceneCmd;
        let t = match self.next() {
            Some(t) => t,
            None => {
                return Err(
                    "parse error: SCENE needs a sub-command (CREATE, CLOSE, TRANSLATE, ROTATE, \
                     ZOOM, HIDE, SHOW, REFRESH, REDRAW, START, STOP, PAUSE, REVERSE, RESET, \
                     SET_TIME_STEP, STATUS, EVENTS)"
                        .into(),
                )
            }
        };
        let cmd = match t.kind {
            TokKind::Keyword(Keyword::Create) => {
                let port = match self.peek() {
                    Some(Token { kind: TokKind::Number(_), .. }) => {
                        let n = self.expect_number("a TCP port")?;
                        if n < 0.0 || n > 65_535.0 || n.fract() != 0.0 {
                            return Err("SCENE CREATE port must be an integer in 0..=65535".into());
                        }
                        n as u16
                    }
                    _ => 0,
                };
                SceneCmd::Create { port }
            }
            TokKind::Keyword(Keyword::Close) => SceneCmd::Close,
            TokKind::Keyword(Keyword::Translate) => {
                self.term(prog)?;
                self.term(prog)?;
                if self.peek().is_some() {
                    self.term(prog)?;
                } else {
                    prog.push(Instr::Push(Value::Num(0.0)));
                }
                SceneCmd::Translate
            }
            TokKind::Keyword(Keyword::Rotate) => {
                self.term(prog)?;
                self.term(prog)?;
                SceneCmd::Rotate
            }
            TokKind::Keyword(Keyword::Zoom) => match self.peek().map(|t| t.kind.clone()) {
                Some(TokKind::Keyword(Keyword::In)) => {
                    self.pos += 1;
                    SceneCmd::ZoomIn
                }
                Some(TokKind::Keyword(Keyword::Out)) => {
                    self.pos += 1;
                    SceneCmd::ZoomOut
                }
                _ => {
                    self.term(prog)?;
                    SceneCmd::Zoom
                }
            },
            TokKind::Keyword(Keyword::Hide) => SceneCmd::Hide(self.scene_which()?),
            TokKind::Keyword(Keyword::Show) => SceneCmd::Show(self.scene_which()?),
            TokKind::Keyword(Keyword::Refresh) => SceneCmd::Refresh,
            TokKind::Keyword(Keyword::Redraw) => SceneCmd::Redraw,
            TokKind::Keyword(Keyword::Start) => SceneCmd::Start,
            TokKind::Keyword(Keyword::Stop) => SceneCmd::Stop,
            TokKind::Keyword(Keyword::Pause) => SceneCmd::Pause,
            TokKind::Keyword(Keyword::Reverse) => SceneCmd::Reverse,
            TokKind::Keyword(Keyword::Reset) => SceneCmd::ResetPlayback,
            TokKind::Keyword(Keyword::SetTimeStep) => {
                self.term(prog)?;
                SceneCmd::SetTimeStep
            }
            TokKind::Keyword(Keyword::Status) => SceneCmd::Status,
            TokKind::Keyword(Keyword::Events) => SceneCmd::Events,
            other => {
                return Err(format!(
                    "parse error at column {}: unknown SCENE sub-command {other} \
                     (expected CREATE, CLOSE, TRANSLATE, ROTATE, ZOOM, HIDE, SHOW, REFRESH, \
                     REDRAW, START, STOP, PAUSE, REVERSE, RESET, SET_TIME_STEP, STATUS or EVENTS)",
                    t.col
                ))
            }
        };
        prog.push(Instr::Scene(cmd));
        Ok(())
    }

    /// `[ NUMBER | "ALL" ]` after HIDE/SHOW — `None` means every object.
    fn scene_which(&mut self) -> Result<Option<usize>, String> {
        match self.peek().map(|t| t.kind.clone()) {
            Some(TokKind::Keyword(Keyword::All)) | None => {
                if self.peek().is_some() {
                    self.pos += 1;
                }
                Ok(None)
            }
            Some(TokKind::Number(_)) => {
                let n = self.expect_number("an object index")?;
                Ok(Some(n as usize))
            }
            Some(other) => Err(format!(
                "parse error: HIDE/SHOW takes an object index or ALL, found {other}"
            )),
        }
    }

    /// `path := IDENT { "." IDENT }` — root `objN` or `system`.
    fn path(&mut self) -> Result<Path, String> {
        let root_name = self.expect_ident("a path root (`objN` or `system`)")?;
        let root = parse_root(&root_name)?;
        self.eat(&TokKind::Dot)?;
        let field = self.expect_field()?;
        let mut comp = None;
        if matches!(self.peek(), Some(t) if t.kind == TokKind::Dot) {
            self.pos += 1;
            let c = self.expect_ident("a component (x, y, z or w)")?;
            comp = Some(match c.to_ascii_lowercase().as_str() {
                "x" => 0usize,
                "y" => 1,
                "z" => 2,
                "w" => 3,
                other => return Err(format!("unknown component `.{other}` (use x, y, z or w)")),
            });
        }
        Ok(Path { root, field, comp })
    }

    fn expr(&mut self, prog: &mut Vec<Instr>) -> Result<(), String> {
        self.term(prog)?;
        loop {
            match self.peek().map(|t| t.kind.clone()) {
                Some(TokKind::Plus) => {
                    self.pos += 1;
                    self.term(prog)?;
                    prog.push(Instr::Add);
                }
                Some(TokKind::Minus) => {
                    self.pos += 1;
                    self.term(prog)?;
                    prog.push(Instr::Sub);
                }
                _ => break,
            }
        }
        Ok(())
    }

    fn term(&mut self, prog: &mut Vec<Instr>) -> Result<(), String> {
        self.unary(prog)?;
        loop {
            match self.peek().map(|t| t.kind.clone()) {
                Some(TokKind::Star) => {
                    self.pos += 1;
                    self.unary(prog)?;
                    prog.push(Instr::Mul);
                }
                Some(TokKind::Slash) => {
                    self.pos += 1;
                    self.unary(prog)?;
                    prog.push(Instr::Div);
                }
                _ => break,
            }
        }
        Ok(())
    }

    fn unary(&mut self, prog: &mut Vec<Instr>) -> Result<(), String> {
        if matches!(self.peek(), Some(t) if t.kind == TokKind::Minus) {
            self.pos += 1;
            self.unary(prog)?;
            prog.push(Instr::Neg);
            return Ok(());
        }
        self.atom(prog)
    }

    fn atom(&mut self, prog: &mut Vec<Instr>) -> Result<(), String> {
        let t = match self.next() {
            Some(t) => t,
            None => return Err("parse error: expected an expression at end of line".into()),
        };
        match t.kind {
            TokKind::Number(n) => {
                prog.push(Instr::Push(Value::Num(n)));
            }
            TokKind::Str(st) => {
                prog.push(Instr::Push(Value::Str(st)));
            }
            TokKind::LBracket => {
                let mut count = 0usize;
                if matches!(self.peek(), Some(t) if t.kind == TokKind::RBracket) {
                    return Err(format!("parse error at column {}: empty vector `[]`", t.col));
                }
                loop {
                    self.expr(prog)?;
                    count += 1;
                    match self.peek().map(|t| t.kind.clone()) {
                        Some(TokKind::Comma) => {
                            self.pos += 1;
                        }
                        _ => break,
                    }
                }
                self.eat(&TokKind::RBracket)?;
                prog.push(Instr::PackList(count));
            }
            TokKind::LParen => {
                self.expr(prog)?;
                self.eat(&TokKind::RParen)?;
            }
            TokKind::Ident(name) => {
                /* builtin call, constant, or path load */
                if matches!(self.peek(), Some(t) if t.kind == TokKind::LParen) {
                    self.pos += 1;
                    let mut argc = 0usize;
                    if !matches!(self.peek(), Some(t) if t.kind == TokKind::RParen) {
                        loop {
                            self.expr(prog)?;
                            argc += 1;
                            match self.peek().map(|t| t.kind.clone()) {
                                Some(TokKind::Comma) => {
                                    self.pos += 1;
                                }
                                _ => break,
                            }
                        }
                    }
                    self.eat(&TokKind::RParen)?;
                    prog.push(Instr::Call(name.to_ascii_lowercase(), argc));
                } else if matches!(self.peek(), Some(t) if t.kind == TokKind::Dot) {
                    /* a dotted path used inside an expression */
                    let root = parse_root(&name)?;
                    self.pos += 1;
                    let field = self.expect_field()?;
                    let mut comp = None;
                    if matches!(self.peek(), Some(t) if t.kind == TokKind::Dot) {
                        self.pos += 1;
                        let c = self.expect_ident("a component (x, y, z or w)")?;
                        comp = Some(match c.to_ascii_lowercase().as_str() {
                            "x" => 0usize,
                            "y" => 1,
                            "z" => 2,
                            "w" => 3,
                            other => {
                                return Err(format!(
                                    "unknown component `.{other}` (use x, y, z or w)"
                                ))
                            }
                        });
                    }
                    prog.push(Instr::Load(Path { root, field, comp }));
                } else {
                    match name.to_ascii_lowercase().as_str() {
                        "pi" => prog.push(Instr::Push(Value::Num(std::f64::consts::PI))),
                        "tau" => prog.push(Instr::Push(Value::Num(std::f64::consts::TAU))),
                        /* a bare identifier: a function parameter or a
                         * LET variable, resolved at execution time */
                        lower => prog.push(Instr::LoadIdent(lower.to_string())),
                    }
                }
            }
            other => {
                return Err(format!(
                    "parse error at column {}: unexpected {other} in expression",
                    t.col
                ));
            }
        }
        Ok(())
    }
}

fn parse_root(name: &str) -> Result<PathRoot, String> {
    let lower = name.to_ascii_lowercase();
    if lower == "system" || lower == "sys" {
        return Ok(PathRoot::System);
    }
    if let Some(idx) = lower.strip_prefix("obj") {
        if let Ok(i) = idx.parse::<usize>() {
            return Ok(PathRoot::Object(i));
        }
    }
    if let Some(idx) = lower.strip_prefix("contact") {
        if let Ok(i) = idx.parse::<usize>() {
            return Ok(PathRoot::Contact(i));
        }
    }
    /* anything else is a USER NAME (an object registered with
     * `NEW ... AS name`, e.g. `dumbell0.m1`), resolved at execution */
    Ok(PathRoot::Named(lower))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_command_compiles_to_postfix() {
        let p = compile_line("set obj0.mass = 2 + 3 * 4").unwrap();
        assert_eq!(
            p,
            vec![
                Instr::Push(Value::Num(2.0)),
                Instr::Push(Value::Num(3.0)),
                Instr::Push(Value::Num(4.0)),
                Instr::Mul,
                Instr::Add,
                Instr::Store(Path {
                    root: PathRoot::Object(0),
                    field: "mass".to_string(),
                    comp: None
                }),
            ]
        );
    }

    #[test]
    fn new_with_inits() {
        let p = compile_line("new sphere { mass = 2, radius = 0.5 }").unwrap();
        assert_eq!(p[0], Instr::NewObject(ShapeKind::Sphere));
        assert_eq!(p[2], Instr::InitField("mass".to_string()));
        assert_eq!(p[4], Instr::InitField("radius".to_string()));
        assert_eq!(p[5], Instr::FinishNew { recompute_inertia: true, name: None });
    }

    #[test]
    fn vector_literal_and_component_path() {
        let p = compile_line("get obj2.position.y").unwrap();
        assert_eq!(
            p,
            vec![Instr::Load(Path {
                root: PathRoot::Object(2),
                field: "position".to_string(),
                comp: Some(1)
            })]
        );
        let p = compile_line("[1, 2, 3]").unwrap();
        assert_eq!(p[3], Instr::PackList(3));
    }

    #[test]
    fn run_and_method() {
        assert_eq!(
            compile_line("run 10 steps 100").unwrap().last().unwrap(),
            &Instr::Run { outputs: 100 }
        );
        assert_eq!(
            compile_line("method sprk leapfrog_2_2 0.001").unwrap()[0],
            Instr::SetMethod(MethodSpec::Sprk {
                table: "ARKODE_SPRK_LEAPFROG_2_2".to_string(),
                dt: 0.001
            })
        );
    }

    #[test]
    fn scene_commands_compile() {
        use crate::vm::SceneCmd;
        assert_eq!(
            compile_line("scene create").unwrap(),
            vec![Instr::Scene(SceneCmd::Create { port: 0 })]
        );
        assert_eq!(
            compile_line("SCENE CREATE 8080").unwrap(),
            vec![Instr::Scene(SceneCmd::Create { port: 8080 })]
        );
        /* translate with an omitted dz gets an implicit 0 */
        let p = compile_line("scene translate 1 2").unwrap();
        assert_eq!(
            p,
            vec![
                Instr::Push(Value::Num(1.0)),
                Instr::Push(Value::Num(2.0)),
                Instr::Push(Value::Num(0.0)),
                Instr::Scene(SceneCmd::Translate),
            ]
        );
        assert_eq!(
            compile_line("scene rotate 15 -5").unwrap().last().unwrap(),
            &Instr::Scene(SceneCmd::Rotate)
        );
        assert_eq!(compile_line("scene zoom in").unwrap(), vec![Instr::Scene(SceneCmd::ZoomIn)]);
        assert_eq!(compile_line("scene zoom out").unwrap(), vec![Instr::Scene(SceneCmd::ZoomOut)]);
        assert_eq!(
            compile_line("scene zoom 2.5").unwrap(),
            vec![Instr::Push(Value::Num(2.5)), Instr::Scene(SceneCmd::Zoom)]
        );
        assert_eq!(compile_line("scene hide all").unwrap(), vec![Instr::Scene(SceneCmd::Hide(None))]);
        assert_eq!(compile_line("scene hide").unwrap(), vec![Instr::Scene(SceneCmd::Hide(None))]);
        assert_eq!(
            compile_line("scene show 2").unwrap(),
            vec![Instr::Scene(SceneCmd::Show(Some(2)))]
        );
        assert_eq!(
            compile_line("scene set_time_step 0.01").unwrap().last().unwrap(),
            &Instr::Scene(SceneCmd::SetTimeStep)
        );
        for (line, cmd) in [
            ("scene start", SceneCmd::Start),
            ("scene stop", SceneCmd::Stop),
            ("scene pause", SceneCmd::Pause),
            ("scene reverse", SceneCmd::Reverse),
            ("scene reset", SceneCmd::ResetPlayback),
            ("scene refresh", SceneCmd::Refresh),
            ("scene redraw", SceneCmd::Redraw),
            ("scene status", SceneCmd::Status),
            ("scene events", SceneCmd::Events),
            ("scene close", SceneCmd::Close),
            ("scene destroy", SceneCmd::Close),
        ] {
            assert_eq!(compile_line(line).unwrap(), vec![Instr::Scene(cmd)], "{line}");
        }
        assert!(compile_line("scene").is_err());
        assert!(compile_line("scene create 99999").is_err());
        assert!(compile_line("scene rotate 15").is_err());
    }

    #[test]
    fn errors_have_positions() {
        let e = compile_line("set = 3").unwrap_err();
        assert!(e.contains("column 5"), "{e}");
        let e = compile_line("get obj0.position.q").unwrap_err();
        assert!(e.contains("component"), "{e}");
        /* a bare identifier now COMPILES (function parameters and LET
         * variables) and resolves at execution instead */
        let p = compile_line("bogusname").unwrap();
        assert_eq!(p, vec![Instr::LoadIdent("bogusname".to_string())]);
    }
}
