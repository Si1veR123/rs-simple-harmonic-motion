#[derive(Debug, Clone)]
pub enum QuadRoot {
    Double(f64, f64),
    Repeated(f64),
    Complex(f64, f64) // first field = a. second field = b. roots are a+bi and a-bi.
}

impl PartialEq for QuadRoot {
    fn eq(&self, other: &Self) -> bool {
        match self {
            QuadRoot::Double(x, y) => {
                if let QuadRoot::Double(a, b) = other {
                    (a == x && b == y) || (a == y && b == x)
                } else {
                    false
                }
            },
            QuadRoot::Repeated(x) => {
                if let QuadRoot::Repeated(y) = other {
                    x == y
                } else {
                    false
                }
            },
            QuadRoot::Complex(x, y) => {
                if let QuadRoot::Complex(a, b) = other {
                    a == x && b.abs() == y.abs()
                } else {
                    false
                }
            },
        }
    }
}

pub fn quadratic_roots(a: f64, b: f64, c: f64) -> QuadRoot {
    let determinant = b.powi(2) - 4.0*a*c;
    
    if determinant < 0.0 {
        let x = -b / (2.0*a);
        let y = (-determinant).sqrt() / (2.0*a);
        QuadRoot::Complex(x, y)
    } else if determinant == 0.0 {
        QuadRoot::Repeated(-b / (2.0*a))
    } else {
        QuadRoot::Double((-b + determinant.sqrt()) / (2.0*a), (-b - determinant.sqrt()) / (2.0*a))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complex_roots() {
        let roots = quadratic_roots(2.0, 4.0, 4.0);
        assert_eq!(roots, QuadRoot::Complex(-1.0, 1.0));
        assert_eq!(roots, QuadRoot::Complex(-1.0, -1.0));
    }

    #[test]
    fn repeated_root() {
        assert_eq!(quadratic_roots(1.0, -8.0, 16.0), QuadRoot::Repeated(4.0))
    }

    #[test]
    fn double_root() {
        let roots = quadratic_roots(2.0, -1.0, -10.0);
        assert_eq!(roots, QuadRoot::Double(2.5, -2.0));
        assert_eq!(roots, QuadRoot::Double(-2.0, 2.5));
    }
}
