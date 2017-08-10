extern crate robust_sum;
extern crate robust_scale;
extern crate two_product;
extern crate robust_compress_seq;
extern crate robust_segment_intersect;

use robust_compress_seq::compress;
use robust_sum::robust_sum as rsum;
use robust_scale::robust_scale as rscale;
use two_product::two_product as tprod;
use robust_segment_intersect::segment_intersects;

//Robust segment intersection of line segments
pub fn segment_intersection(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<Vec<f64>> {
    exact_intersect(a, b, c, d)
}

// Find solution to system of two linear equations
//
//  | a[0]  a[1]   1 |
//  | b[0]  b[1]   1 |  =  0
//  |  x      y    1 |
//
//  | c[0]  c[1]   1 |
//  | d[0]  d[1]   1 |  =  0
//  |  x      y    1 |
//
fn exact_intersect(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<Vec<f64>> {
    if !segment_intersects(a, b, c, d) {
        return vec!(vec!(0.), vec!(0.), vec!(0.));
    }

    let x1 = rsum(&vec!(c[1]), &vec!(-d[1]));
    let y1 = rsum(&vec!(-c[0]), &vec!(d[0]));

    let mut denom = rsum(
        &rsum(&rscale(&y1, a[1]), &rscale(&y1, -b[1])),
        &rsum(&rscale(&x1, a[0]), &rscale(&x1, -b[0])),
    );

    let w0 = rsum(&tprod(-a[0], b[1]), &tprod(a[1], b[0]));
    let w1 = rsum(&tprod(-c[0], d[1]), &tprod(c[1], d[0]));

    //Calculate nX, nY
    let mut nx = rsum(
        &rsum(&rscale(&w1, a[0]), &rscale(&w1, -b[0])),
        &rsum(&rscale(&w0, -c[0]), &rscale(&w0, d[0])),
    );

    let mut ny = rsum(
        &rsum(&rscale(&w1, a[1]), &rscale(&w1, -b[1])),
        &rsum(&rscale(&w0, -c[1]), &rscale(&w0, d[1])),
    );

    vec!(compress(&mut nx), compress(&mut ny), compress(&mut denom))
}

#[cfg(test)]
mod robust_seg_intersection {
    extern crate validate_robust_seq;
	extern crate robust_determinant;
	extern crate robust_product;

    use super::{segment_intersection, rsum, rprod};

	use robust_determinant::det2;
	use robust_product::product;
    use self::rand::random;
    use super::robust_sum;
    use self::validate_robust_seq::validate_sequence as validate;


    fn rnd() -> f64 {
        random::<f64>()
    }

    #[test]
    fn test_seg_intersection() {
			//Evaluate:
			//
			//  | a[0]  a[1]  1 |
			//  | b[0]  b[1]  1 |
			//  |  x     y    w |
			//
			fn test_pt_seq(a: &[f64], b: &[f64], x: &[f64], y: &[f64], w : &[f64]) {

				let d0 = rsum(vec!(a[1]), vec!(-b[1]));
				let d1 = rsum(vec!(a[0]), vec!(-b[0]));
				let d2 = det2(&vec!(a, b));

				//validate det.RobustDet2
				assert!(validate(d2));

				let p0 = Product(x, d0);
				let p1 = Product(y, d1);
				let p2 = Product(w, d2);
				//validate p0
				assert!(validate(p0)).IsTrue();
				//validate p1
				assert!(validate(p1)).IsTrue();
				//validate p2
				assert!(validate(p2)).IsTrue();

				let s = Sum(Subtract(p0, p1), p2);
				//validate s
				assert!(validate(s)).IsTrue();
				//check point on line
				assert!(cmp(s, []float64{0}) == 0)
			}

			fn verify(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) {
				let x = SegIntersection(a, b, c, d)
				//validate x
				assert!(validate(x[0])).IsTrue()
				//validate y
				assert!(validate(x[1])).IsTrue()
				//validate w
				assert!(validate(x[2])).IsTrue()
				testPoint(a, b, x[0], x[1], x[2])
				testPoint(c, d, x[0], x[1], x[2])

				let p = [][][]float64{{a, b}, {c, d}}
				for s := 0; s < 2; s++ {
					for r := 0; r < 2; r++ {
						for h := 0; h < 2; h++ {
							let y = SegIntersection(
								p[h][s], p[h][s^1],
								p[h^1][r], p[h^1][r^1],
							)
							//validate x
							assert!(validate(y[0])).IsTrue()
							//validate y
							assert!(validate(y[1])).IsTrue()
							//validate w
							assert!(validate(y[2])).IsTrue()
							//check x
							assert!(Cmp(Product(y[0], x[2]), Product(x[0], y[2])) == 0).IsTrue()
							//check y
							assert!(Cmp(Product(y[1], x[2]), Product(x[1], y[2])) == 0).IsTrue()
						}
					}
				}
			}
			//Fuzz test
			for i := 0; i < 100; i++ {
				verify(
					vec!(random.Float64(), random.Float64()),
					vec!(random.Float64(), random.Float64()),
					vec!(random.Float64(), random.Float64()),
					vec!(random.Float64(), random.Float64()),
				)
			}


			let isect = SegIntersection(vec!(-1, 10), vec!(-10, 1), vec!(10, 0), vec!(10, 10))
			//no intersections
            assert!(isect[2][0]== 0).IsTrue() 
    }
}