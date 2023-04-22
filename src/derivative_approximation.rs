
// used test to find a rough optimal
const DEFAULT_DELTA: f64 = 0.00045;

pub fn simple_gradient_approx_delta<F>(f: F, x: f64, delta: f64) -> f64
    where F: Fn(f64) -> f64
{
    let pre_y = f(x - delta);
    let post_y = f(x + delta);

    let gradient = (post_y - pre_y) / (2.0*delta);
    gradient
}

pub fn simple_gradient_approx<F>(f: F, x: f64) -> f64
    where F: Fn(f64) -> f64
{
    simple_gradient_approx_delta(f, x, DEFAULT_DELTA)
}

#[cfg(test)]
mod tests {
    use super::*;

    // within this percentage
    const MAX_ERROR: f64 = 0.0000001;

    fn run_test_function(f: fn(f64) -> f64, d: fn(f64) -> f64) {
        for x in -10000..=10000 {
            let calculated_gradient = simple_gradient_approx(f, x as f64);
            let actual_gradient = d(x as f64);

            // if near 0 (likely 0 with rounding error), dont use % error calculation to avoid divide by 0.
            if actual_gradient.abs() < 0.000001  {
                if calculated_gradient.abs() > 0.000001 {
                    panic!("Panicked at gradient 0. Calculated {}.", calculated_gradient)
                }
                // else accurate approximation

            } else {
                let percent_error = (actual_gradient-calculated_gradient)/actual_gradient;
                
                if percent_error.abs() > MAX_ERROR {
                    panic!("ERROR is {}% for x = {}. Calculated {}. Actual {}.\n", percent_error.abs(), x, calculated_gradient, actual_gradient);
                }
            }
        }
    }

    #[test]
    fn quadratic() {
        run_test_function(|x| x.powi(2), |x| x*2.0);
    }

    #[test]
    fn sin() {
        run_test_function(f64::sin, f64::cos);
    }

    #[test]
    fn flat() {
        for x in -1000..1000 {
            let calc = simple_gradient_approx(|_x| 5.0, x as f64);
            let actual = 0.0;
            assert_eq!(calc, actual);
        }
    }

    #[test]
    fn exp() {
        run_test_function(f64::exp, f64::exp)
    }

    #[test]
    fn cubic() {
        run_test_function(|x| x.powi(3), |x| 3.0*(x.powi(2)))
    }

    fn find_min_error() -> f64 {
        let funcs: Vec<(fn(f64) -> f64, fn(f64) -> f64)> = vec![
            (|x: f64| x.powi(2), |x| x*2.0),
            (|x: f64| x.powi(4), |x| x*4.0),
            (|x: f64| x.powi(10), |x| x*10.0),
            (|x: f64| x.sin(), |x: f64| x.cos()),
            (|x: f64| x.cos(), |x: f64| -x.sin()),
            (|x: f64| x.exp(), |x: f64| x.exp()),
            (|x: f64| 2.0f64.powf(x), |x: f64| 2.0f64.ln() * 2.0f64.powf(x)),
            (|_x: f64| 5.0, |_x: f64| 0.0),
        ];

        let mut delta = 10.0;
        let delta_adjustment = 0.999;
        let n = 10000;
        let mut best_error = f64::INFINITY;
        let mut best_delta = 0.0;

        for _ in 0..n {
            let mut error = 0.0;
            for f in funcs.iter() {
                for x in -500..=500 {
                    let calculated_gradient = simple_gradient_approx_delta(f.0, x as f64, delta);
                    let actual_gradient = (f.1)(x as f64);
                    let gradient_error = actual_gradient - calculated_gradient;

                    error += gradient_error.abs();
                }
            }

            if error < best_error {
                best_error = error;
                best_delta = delta;
            }

            delta *= delta_adjustment;
        }

        best_delta
    }

    #[ignore]
    #[test]
    fn find_min_error_repeated() {
        let mut delta = 0.0;
        for _ in 0..10 {
            delta += find_min_error()
        }
        println!("{}", delta/10.0)
    }
}
