use crate::quadratic_roots::{quadratic_roots, QuadRoot};

pub enum Return {
    Displacement,
    Velocity,
    Acceleration
}

pub enum InitialConditions {
    DisplacementVelocity(f64, f64),
    DisplacementAcceleration(f64, f64),
    VelocityAcceleration(f64, f64)
}

pub fn second_order_homogenous_solver(a: f64, b: f64, c: f64, t: f64, initial: InitialConditions, return_value: Return) -> f64 {
    // solve a differential in the form ax'' + bx' + cx = 0 at t, given initial conditions
    // returns displacement, velocity or acceleration at t depending on return_value

    let roots = quadratic_roots(a, b, c);

    let constant_a: f64;
    let constant_b: f64;
    match roots {
        QuadRoot::Double(root_a, root_b) => {
            match initial {
                // simultaneous equations when t=0, to find A and B
                InitialConditions::DisplacementVelocity(x, v) => {
                    constant_b = (v - (root_a*x)) / (root_b - root_a);
                    constant_a = x - constant_b;
                },
                InitialConditions::DisplacementAcceleration(x, a) => {
                    constant_b = (a - root_a.powi(2)*x) / (root_b.powi(2) - root_a.powi(2));
                    constant_a = x - constant_b;
                },
                InitialConditions::VelocityAcceleration(v, a) => {
                    constant_b = (a - root_a*v) / (root_b.powi(2) - root_b*root_a);
                    constant_a = (v - root_b*constant_b) / root_a;
                }
            }

            match return_value {
                Return::Displacement => return constant_a*(root_a*t).exp() + constant_b*(root_b*t).exp(),
                Return::Velocity => return root_a*constant_a*(root_a*t).exp() + root_b*constant_b*(root_b*t).exp(),
                Return::Acceleration => return root_a.powi(2)*constant_a*(root_a*t).exp() + root_b.powi(2)*constant_b*(root_b*t).exp(),
            }
        }
        QuadRoot::Repeated(root) => {
            match initial {
                // simultaneous equations when t=0, to find A and B
                InitialConditions::DisplacementVelocity(x, v) => {
                    constant_a = x;
                    constant_b = v - (root * x);
                },
                InitialConditions::DisplacementAcceleration(x, a) => {
                    constant_a = x;
                    constant_b = (a - (root.powi(2)*x)) / (2.0 * root);
                },
                InitialConditions::VelocityAcceleration(v, a) => {
                    constant_b = (a - root*v) / root;
                    constant_a = (v - constant_b) / root;
                }
            }

            match return_value {
                Return::Displacement => return (root*t).exp() * (constant_a + constant_b*t),
                Return::Velocity => return (root*t).exp() * (root * (constant_a + constant_b*t) + constant_b),
                Return::Acceleration => return (root*t).exp() * (root.powi(2)*(constant_a+t*constant_b) + 2.0*constant_b*root)
            }
        },
        QuadRoot::Complex(p, q) => {
            match initial {
                InitialConditions::DisplacementVelocity(x, v) => {
                    constant_a = x;
                    constant_b = (v - (x*p)) / q;
                },
                InitialConditions::DisplacementAcceleration(x, a) => {
                    constant_a = x;
                    constant_b = (a - x*(p.powi(2)-q.powi(2))) / (2.0*p*q);
                },
                InitialConditions::VelocityAcceleration(v, a) => {
                    constant_a = - (a - 2.0*p*v) / (p.powi(2) + q.powi(2));
                    constant_b = (v - constant_a*p) / q;
                },
            }
            
            match return_value {
                Return::Displacement => return (p*t).exp() * (constant_a*(q*t).cos() + constant_b*(q*t).sin()),
                Return::Velocity => return (p*t).exp() * ((constant_a*p + constant_b*q)*(q*t).cos() + (constant_b*p - constant_a*q)*(q*t).sin()),
                Return::Acceleration => {
                    let sq_diff = p.powi(2) - q.powi(2);
                    return (p*t).exp() * ( (q*t).cos() * (constant_a*sq_diff + 2.0*constant_b*q*p) + (q*t).sin() * (constant_b*sq_diff - 2.0*constant_a*q*p) )
                }
            }
        },
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_roots() {
        {
            // initial conditions: displacement and velocity
            let displacement = second_order_homogenous_solver(1.0, 5.0, 6.0, 5.0, InitialConditions::DisplacementVelocity(2.0, 4.0), Return::Displacement);
            let error_d = (displacement - 0.00045155).abs();
            assert!(error_d < 0.00000001);

            let velocity = second_order_homogenous_solver(1.0, 5.0, 6.0, 5.0, InitialConditions::DisplacementVelocity(2.0, 4.0), Return::Velocity);
            let error_v = (velocity - -0.00090065).abs();
            assert!(error_v < 0.00000001);

            let acceleration  = second_order_homogenous_solver(1.0, 5.0, 6.0, 5.0, InitialConditions::DisplacementVelocity(2.0, 4.0), Return::Acceleration);
            let error_a = (acceleration - 0.00179397).abs();
            assert!(error_a < 0.00000001);
        }

        {
            // initial conditions: displacement and acceleration
            let displacement = second_order_homogenous_solver(1.0, 5.0, 6.0, 2.0, InitialConditions::DisplacementAcceleration(2.0, 4.0), Return::Displacement);
            let error_d = (displacement - 0.04930079).abs();
            assert!(error_d < 0.00000001);

            let velocity = second_order_homogenous_solver(1.0, 5.0, 6.0, 2.0, InitialConditions::DisplacementAcceleration(2.0, 4.0), Return::Velocity);
            let error_v = (velocity - -0.09661857).abs();
            assert!(error_v < 0.00000001);

            let acceleration = second_order_homogenous_solver(1.0, 5.0, 6.0, 2.0, InitialConditions::DisplacementAcceleration(2.0, 4.0), Return::Acceleration);
            let error_a = (acceleration - 0.18728813).abs();
            assert!(error_a < 0.00000001);
        }

        {
            // initial conditions: velocity and acceleration
            let displacement = second_order_homogenous_solver(1.0, 5.0, 6.0, 3.0, InitialConditions::VelocityAcceleration(2.0, 4.0), Return::Displacement);
            println!("{displacement}");
            let error_d = (displacement - -0.01206467).abs();
            assert!(error_d < 0.00000001);

            let velocity = second_order_homogenous_solver(1.0, 5.0, 6.0, 3.0, InitialConditions::VelocityAcceleration(2.0, 4.0), Return::Velocity);
            let error_v = (velocity - 0.02380024).abs();
            assert!(error_v < 0.00000001);

            let acceleration = second_order_homogenous_solver(1.0, 5.0, 6.0, 3.0, InitialConditions::VelocityAcceleration(2.0, 4.0), Return::Acceleration);
            let error_a = (acceleration - -0.04661321).abs();
            assert!(error_a < 0.00000001);
        }
    }

    #[test]
    fn single_root() {
        {
            // initial conditions: displacement and velocity
            let displacement = second_order_homogenous_solver(1.0, -8.0, 16.0, 2.0, InitialConditions::DisplacementVelocity(2.0, 4.0), Return::Displacement);
            let error_d = (displacement - -17885.74792).abs();
            assert!(error_d < 0.00001);

            let velocity = second_order_homogenous_solver(1.0, -8.0, 16.0, 2.0, InitialConditions::DisplacementVelocity(2.0, 4.0), Return::Velocity);
            let error_v = (velocity - -83466.82364).abs();
            assert!(error_v < 0.00001);

            let acceleration = second_order_homogenous_solver(1.0, -8.0, 16.0, 2.0, InitialConditions::DisplacementVelocity(2.0, 4.0), Return::Acceleration);
            let error_a = (acceleration - -381562.6223).abs();
            assert!(error_a < 0.0001);
        }

        {
            // initial conditions: displacement and acceleration
            let displacement = second_order_homogenous_solver(1.0, -8.0, 16.0, 2.0, InitialConditions::DisplacementAcceleration(2.0, 4.0), Return::Displacement);
            let error_d = (displacement - -14904.78994).abs();
            assert!(error_d < 0.00001);

            let velocity = second_order_homogenous_solver(1.0, -8.0, 16.0, 2.0, InitialConditions::DisplacementAcceleration(2.0, 4.0), Return::Velocity);
            let error_v = (velocity - -70052.5127).abs();
            assert!(error_v < 0.0001);

            let acceleration = second_order_homogenous_solver(1.0, -8.0, 16.0, 2.0, InitialConditions::DisplacementAcceleration(2.0, 4.0), Return::Acceleration);
            let error_a = (acceleration - -321943.4626).abs();
            assert!(error_a < 0.0001);
        }

        {
            // initial conditions: velocity and acceleration
            let displacement = second_order_homogenous_solver(1.0, -8.0, 16.0, 2.0, InitialConditions::VelocityAcceleration(2.0, 4.0), Return::Displacement);
            let error_d = (displacement - -3726.197484).abs();
            assert!(error_d < 0.000001);

            let velocity = second_order_homogenous_solver(1.0, -8.0, 16.0, 2.0, InitialConditions::VelocityAcceleration(2.0, 4.0), Return::Velocity);
            let error_v = (velocity - -17885.74792).abs();
            assert!(error_v < 0.00001);

            let acceleration = second_order_homogenous_solver(1.0, -8.0, 16.0, 2.0, InitialConditions::VelocityAcceleration(2.0, 4.0), Return::Acceleration);
            let error_a = (acceleration - -83466.82364).abs();
            assert!(error_a < 0.00001);
        }
    }

    #[test]
    fn no_roots() {
        {
            // initial conditions: displacement and velocity
            let displacement = second_order_homogenous_solver(2.0, 4.0, 4.0, 2.0, InitialConditions::DisplacementVelocity(2.0, 4.0), Return::Displacement);
            let error_d = (displacement - 0.6257214489).abs();
            assert!(error_d < 0.0000000001);

            let velocity = second_order_homogenous_solver(2.0, 4.0, 4.0, 2.0, InitialConditions::DisplacementVelocity(2.0, 4.0), Return::Velocity);
            let error_v = (velocity - -1.209757598).abs();
            assert!(error_v < 0.000000001);

            let acceleration = second_order_homogenous_solver(2.0, 4.0, 4.0, 2.0, InitialConditions::DisplacementVelocity(2.0, 4.0), Return::Acceleration);
            println!("{acceleration}");
            let error_a = (acceleration - 1.168072299).abs();
            assert!(error_a < 0.000000001);
        }

        {
            // initial conditions: displacement and acceleration
            let displacement = second_order_homogenous_solver(2.0, 4.0, 4.0, 2.0, InitialConditions::DisplacementAcceleration(2.0, 4.0), Return::Displacement);
            let error_d = (displacement - -0.3587587496).abs();
            assert!(error_d < 0.0000000001);

            let velocity = second_order_homogenous_solver(2.0, 4.0, 4.0, 2.0, InitialConditions::DisplacementAcceleration(2.0, 4.0), Return::Velocity);
            let error_v = (velocity - 0.2252774).abs();
            assert!(error_v < 0.0000001);

            let acceleration = second_order_homogenous_solver(2.0, 4.0, 4.0, 2.0, InitialConditions::DisplacementAcceleration(2.0, 4.0), Return::Acceleration);
            let error_a = (acceleration - 0.2669626993).abs();
            assert!(error_a < 0.0000000001);
        }

        {
            // initial conditions: velocity and acceleration
            let displacement = second_order_homogenous_solver(2.0, 4.0, 4.0, 2.0, InitialConditions::VelocityAcceleration(2.0, 4.0), Return::Displacement);
            let error_d = (displacement - -0.02084264964).abs();
            assert!(error_d < 0.00000000001);

            let velocity = second_order_homogenous_solver(2.0, 4.0, 4.0, 2.0, InitialConditions::VelocityAcceleration(2.0, 4.0), Return::Velocity);
            let error_v = (velocity - 0.6257214489).abs();
            assert!(error_v < 0.0000000001);

            let acceleration = second_order_homogenous_solver(2.0, 4.0, 4.0, 2.0, InitialConditions::VelocityAcceleration(2.0, 4.0), Return::Acceleration);
            let error_a = (acceleration - -1.209757598).abs();
            assert!(error_a < 0.000000001);
        }
    }
}
