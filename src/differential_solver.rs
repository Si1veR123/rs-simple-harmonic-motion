use crate::quadratic_roots::{quadratic_roots, QuadRoot};


pub fn second_order_homogenous_solver(a: f64, b: f64, c: f64, t: f64) {
    let roots = quadratic_roots(a, b, c);
    match roots {
        QuadRoot::Double(x, y) => {
            todo!()
        }
        QuadRoot::Repeated(x) => todo!(),
        QuadRoot::Complex(x, y) => todo!(),
    }
}
