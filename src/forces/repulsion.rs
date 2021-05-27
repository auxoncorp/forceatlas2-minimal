use crate::{iter::*, layout::*, util::*};

use itertools::izip;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use num_traits::{One, Zero};
use std::ops::{SubAssign, AddAssign};

pub fn apply_repulsion(layout: &mut Layout) {
	let kr = layout.settings.kr.clone();
	let mut di = valloc(layout.settings.dimensions);
	for Node {
		mass: n1_mass,
		n2_iter,
		pos: n1_pos,
		speed: n1_speed,
		..
	} in layout.iter_nodes()
	{
		let n1_mass = n1_mass.clone() + Coord::one();
		for Node2 {
			mass: n2_mass,
			pos: n2_pos,
			speed: n2_speed,
			..
		} in n2_iter
		{
			di.clone_from_slice(n2_pos);

			let d2 = di
				.iter_mut()
				.zip(n1_pos.iter())
				.map(|(di, n1_pos)| {
					*di -= n1_pos.clone();
					(*di).powi(2)
				})
				.sum::<Coord>();
			if d2.is_zero() {
				continue;
			}

			let f = n1_mass.clone() * (n2_mass.clone() + Coord::one()) / d2 * kr.clone();

			izip!(n1_speed.iter_mut(), n2_speed.iter_mut(), di.iter()).for_each(
				|(n1_speed, n2_speed, di)| {
					let s = f.clone() * di.clone();
					*n1_speed -= s.clone();
					*n2_speed += s;
				},
			);
		}
	}
}

#[cfg(feature = "parallel")]
pub fn apply_repulsion_parallel(layout: &mut Layout) {
	let kr = layout.settings.kr.clone();
	let dimensions = layout.settings.dimensions;

	for chunk_iter in layout.iter_par_nodes(layout.settings.chunk_size.unwrap()) {
		chunk_iter.for_each(|n1_iter| {
			let mut di = valloc(dimensions);
			for n1 in n1_iter {
				let n1_mass = n1.mass.clone() + Coord::one();
				for n2 in n1.n2_iter {
					di.clone_from_slice(n2.pos);

					let d2 = di
						.iter_mut()
						.zip(n1.pos.iter())
						.map(|(di, n1_pos)| {
							*di -= n1_pos.clone();
							(*di).powi(2)
						})
						.sum::<Coord>();
					if d2.is_zero() {
						continue;
					}

					let f = n1_mass.clone() * (n2.mass.clone() + Coord::one()) / d2 * kr.clone();

					izip!(n1.speed.iter_mut(), n2.speed.iter_mut(), di.iter()).for_each(
						|(n1_speed, n2_speed, di)| {
							let s = f.clone() * di.clone();
							*n1_speed -= s.clone();
							*n2_speed += s;
						},
					);
				}
			}
		});
	}
}

pub fn apply_repulsion_2d(layout: &mut Layout) {
	let kr = layout.settings.kr;
	for Node {
		mass: n1_mass,
		n2_iter,
		pos: n1_pos,
		speed: n1_speed,
		..
	} in layout.iter_nodes()
	{
		let n1_mass = *n1_mass + Coord::one();
		for Node2 {
			mass: n2_mass,
			pos: n2_pos,
			speed: n2_speed,
			..
		} in n2_iter
		{
			let dx = unsafe { *n2_pos.get_unchecked(0) - *n1_pos.get_unchecked(0) };
			let dy = unsafe { *n2_pos.get_unchecked(1) - *n1_pos.get_unchecked(1) };

			let d2 = dx * dx + dy * dy;
			if d2.is_zero() {
				continue;
			}

			let f = n1_mass * (*n2_mass + Coord::one()) / d2 * kr;

			let vx = f * dx;
			let vy = f * dy;
			unsafe { n1_speed.get_unchecked_mut(0) }.sub_assign(vx); // n1_speed[0] -= f * dx
			unsafe { n1_speed.get_unchecked_mut(1) }.sub_assign(vy); // n1_speed[1] -= f * dy
			unsafe { n2_speed.get_unchecked_mut(0) }.add_assign(vx); // n2_speed[0] += f * dx
			unsafe { n2_speed.get_unchecked_mut(1) }.add_assign(vy); // n2_speed[1] += f * dy
		}
	}
}

#[cfg(feature = "parallel")]
pub fn apply_repulsion_2d_parallel(
	layout: &mut Layout,
) {
	let kr = layout.settings.kr;
	for chunk_iter in layout.iter_par_nodes(layout.settings.chunk_size.unwrap()) {
		chunk_iter.for_each(|n1_iter| {
			for n1 in n1_iter {
				let n1_mass = *n1.mass + Coord::one();
				for n2 in n1.n2_iter {
					let dx = unsafe { *n2.pos.get_unchecked(0) - *n1.pos.get_unchecked(0) };
					let dy = unsafe { *n2.pos.get_unchecked(1) - *n1.pos.get_unchecked(1) };

					let d2 = dx * dx + dy * dy;
					if d2.is_zero() {
						continue;
					}

					let f = n1_mass * (*n2.mass + Coord::one()) / d2 * kr;

					let vx = f * dx;
					let vy = f * dy;
					unsafe { n1.speed.get_unchecked_mut(0) }.sub_assign(vx); // n1_speed[0] -= f * dx
					unsafe { n1.speed.get_unchecked_mut(1) }.sub_assign(vy); // n1_speed[1] -= f * dy
					unsafe { n2.speed.get_unchecked_mut(0) }.add_assign(vx); // n2_speed[0] += f * dx
					unsafe { n2.speed.get_unchecked_mut(1) }.add_assign(vy); // n2_speed[1] += f * dy
				}
			}
		});
	}
}


#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub fn apply_repulsion_2d_simd_f32(layout: &mut Layout) {
	#[cfg(target_arch = "x86")]
	use std::arch::x86::*;
	#[cfg(target_arch = "x86_64")]
	use std::arch::x86_64::*;

	for (n1, (n1_mass, n1_pos_s)) in layout.masses.iter().zip(layout.points.iter()).enumerate() {
		let mut n2_iter = layout.points.iter();
		let n1_mass = n1_mass + 1.0f32;
		let n1_speed = unsafe {
			layout
				.speeds
				.points
				.as_mut_ptr()
				.add(n1 * layout.settings.dimensions)
		};
		let n1_pos = unsafe {
			_mm256_set_ps(
				*n1_pos_s.get_unchecked(1),
				*n1_pos_s.get_unchecked(0),
				*n1_pos_s.get_unchecked(1),
				*n1_pos_s.get_unchecked(0),
				*n1_pos_s.get_unchecked(1),
				*n1_pos_s.get_unchecked(0),
				*n1_pos_s.get_unchecked(1),
				*n1_pos_s.get_unchecked(0),
			)
		};

		/*unsafe {
			assert_eq!(std::mem::transmute::<__m256, (f32,f32,f32,f32,f32,f32,f32,f32)>(n1_pos), (
				layout.points.get(n1)[0],
				layout.points.get(n1)[1],
				layout.points.get(n1)[0],
				layout.points.get(n1)[1],
				layout.points.get(n1)[0],
				layout.points.get(n1)[1],
				layout.points.get(n1)[0],
				layout.points.get(n1)[1],
			));
		}*/

		// This loop iterates on nodes by 4
		let n2_max = n1 & (usize::MAX - 3);
		let mut n2 = 0usize;
		while n2 < n2_max {
			unsafe {
				// [n2_x, n2_y, n3_x, n3_y, n4_x, n4_y, n5_x, n5_y]
				let n2345_pos = _mm256_loadu_ps(n2_iter.next_unchecked(4));

				/*assert_eq!(std::mem::transmute::<__m256, (f32,f32,f32,f32,f32,f32,f32,f32)>(n2345_pos), (
					layout.points.get(n2)[0],
					layout.points.get(n2)[1],
					layout.points.get(n2+1)[0],
					layout.points.get(n2+1)[1],
					layout.points.get(n2+2)[0],
					layout.points.get(n2+2)[1],
					layout.points.get(n2+3)[0],
					layout.points.get(n2+3)[1],
				));*/

				// [dx(n1,n2), dy(n1,n2), dx(n1,n3), dy(n1,n3), dx(n1,n4), dy(n1,n4), dx(n1,n5), dy(n1,n5)]
				let dxy = _mm256_sub_ps(n2345_pos, n1_pos);

				/*assert_eq!(std::mem::transmute::<__m256, (f32,f32,f32,f32,f32,f32,f32,f32)>(dxy), (
					layout.points.get(n2)[0] - layout.points.get(n1)[0],
					layout.points.get(n2)[1] - layout.points.get(n1)[1],
					layout.points.get(n2+1)[0] - layout.points.get(n1)[0],
					layout.points.get(n2+1)[1] - layout.points.get(n1)[1],
					layout.points.get(n2+2)[0] - layout.points.get(n1)[0],
					layout.points.get(n2+2)[1] - layout.points.get(n1)[1],
					layout.points.get(n2+3)[0] - layout.points.get(n1)[0],
					layout.points.get(n2+3)[1] - layout.points.get(n1)[1],
				));*/

				// ([dx(n1,n2)^2, dx(n1,n3)^2, dx(n1,n4)^2, dx(n1,n5)^2], [dy(n1,n2)^2, dy(n1,n3)^2, dy(n1,n4)^2, dy(n1,n5)^2])
				let (dx2, dy2): (__m128, __m128) = std::mem::transmute(_mm256_permutevar8x32_ps(
					_mm256_mul_ps(dxy, dxy),
					_mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0),
				));

				/*assert_eq!((std::mem::transmute::<__m128, (f32,f32,f32,f32)>(dx2), std::mem::transmute::<__m128, (f32,f32,f32,f32)>(dy2)), (
					(
						(layout.points.get(n2)[0] - layout.points.get(n1)[0]).powi(2),
						(layout.points.get(n2+1)[0] - layout.points.get(n1)[0]).powi(2),
						(layout.points.get(n2+2)[0] - layout.points.get(n1)[0]).powi(2),
						(layout.points.get(n2+3)[0] - layout.points.get(n1)[0]).powi(2)
					),
					(
						(layout.points.get(n2)[1] - layout.points.get(n1)[1]).powi(2),
						(layout.points.get(n2+1)[1] - layout.points.get(n1)[1]).powi(2),
						(layout.points.get(n2+2)[1] - layout.points.get(n1)[1]).powi(2),
						(layout.points.get(n2+3)[1] - layout.points.get(n1)[1]).powi(2)
					)
				));*/

				// [d(n1,n2), d(n1,n3), d(n1,n4), d(n1,n5)]
				let d2 = _mm_add_ps(dx2, dy2);
				// TODO maybe check zero

				/*assert_eq!(std::mem::transmute::<__m128, (f32,f32,f32,f32)>(d2), (
					(layout.points.get(n2)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2)[1] - layout.points.get(n1)[1]).powi(2),
					(layout.points.get(n2+1)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+1)[1] - layout.points.get(n1)[1]).powi(2),
					(layout.points.get(n2+2)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+2)[1] - layout.points.get(n1)[1]).powi(2),
					(layout.points.get(n2+3)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+3)[1] - layout.points.get(n1)[1]).powi(2),
				));*/

				let degs2 = n1_mass * (layout.masses.get_unchecked(n2) + 1.);
				let degs3 = n1_mass * (layout.masses.get_unchecked(n2 + 1) + 1.);
				let degs4 = n1_mass * (layout.masses.get_unchecked(n2 + 2) + 1.);
				let degs5 = n1_mass * (layout.masses.get_unchecked(n2 + 3) + 1.);
				let f = _mm_mul_ps(
					_mm_div_ps(_mm_set_ps(degs5, degs4, degs3, degs2), d2),
					_mm_set_ps(
						layout.settings.kr,
						layout.settings.kr,
						layout.settings.kr,
						layout.settings.kr,
					),
				);

				/*assert_eq!(std::mem::transmute::<__m128, (f32,f32,f32,f32)>(f), (
					degs2/((layout.points.get(n2)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr,
					degs3/((layout.points.get(n2+1)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+1)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr,
					degs4/((layout.points.get(n2+2)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+2)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr,
					degs5/((layout.points.get(n2+3)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+3)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr,
				));*/

				let n2345_speed = layout
					.speeds
					.points
					.as_mut_ptr()
					.add(n2 * layout.settings.dimensions);
				let (n1_speed_v, n2345_speed_v): (__m128, __m256) =
					(_mm_loadu_ps(n1_speed), _mm256_loadu_ps(n2345_speed));

				let fxy = _mm256_mul_ps(
					dxy,
					_mm256_permutevar8x32_ps(
						_mm256_set_m128(f, f),
						_mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0),
					),
				);

				/*assert_eq!(std::mem::transmute::<__m256, (f32,f32,f32,f32,f32,f32,f32,f32)>(fxy), (
					degs2/((layout.points.get(n2)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr*(layout.points.get(n2)[0] - layout.points.get(n1)[0]),
					degs2/((layout.points.get(n2)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr*(layout.points.get(n2)[1] - layout.points.get(n1)[1]),
					degs3/((layout.points.get(n2+1)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+1)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr*(layout.points.get(n2+1)[0] - layout.points.get(n1)[0]),
					degs3/((layout.points.get(n2+1)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+1)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr*(layout.points.get(n2+1)[1] - layout.points.get(n1)[1]),
					degs4/((layout.points.get(n2+2)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+2)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr*(layout.points.get(n2+2)[0] - layout.points.get(n1)[0]),
					degs4/((layout.points.get(n2+2)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+2)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr*(layout.points.get(n2+2)[1] - layout.points.get(n1)[1]),
					degs5/((layout.points.get(n2+3)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+3)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr*(layout.points.get(n2+3)[0] - layout.points.get(n1)[0]),
					degs5/((layout.points.get(n2+3)[0] - layout.points.get(n1)[0]).powi(2)+
					(layout.points.get(n2+3)[1] - layout.points.get(n1)[1]).powi(2))*layout.settings.kr*(layout.points.get(n2+3)[1] - layout.points.get(n1)[1]),
				));*/

				_mm256_storeu_ps(n2345_speed, _mm256_add_ps(n2345_speed_v, fxy));
				_mm_storeu_ps(
					n1_speed,
					_mm_sub_ps(
						n1_speed_v,
						_mm_add_ps(_mm256_extractf128_ps(fxy, 1), _mm256_extractf128_ps(fxy, 0)),
					),
				);
			}

			n2 += 4;
		}

		// Remaining iterations (if n1 is not multiple of 4)
		while n2 < n1 {
			let n2_pos = unsafe { layout.points.get_unchecked(n2) };

			let dx = unsafe { *n2_pos.get_unchecked(0) - *n1_pos_s.get_unchecked(0) };
			let dy = unsafe { *n2_pos.get_unchecked(1) - *n1_pos_s.get_unchecked(1) };

			let d2 = dx * dx + dy * dy;
			if d2.is_zero() {
				continue;
			}

			let f = (n1_mass * (unsafe { layout.masses.get_unchecked(n2) } + 1.)) / d2
				* layout.settings.kr;

			let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);
			let vx = f * dx;
			let vy = f * dy;
			unsafe { n1_speed.get_unchecked_mut(0) }.sub_assign(vx); // n1_speed[0] -= f * dx
			unsafe { n1_speed.get_unchecked_mut(1) }.sub_assign(vy); // n1_speed[1] -= f * dy
			unsafe { n2_speed.get_unchecked_mut(0) }.add_assign(vx); // n2_speed[0] += f * dx
			unsafe { n2_speed.get_unchecked_mut(1) }.add_assign(vy); // n2_speed[1] += f * dy

			n2 += 1;
		}
	}
}


#[cfg(all(feature = "parallel", any(target_arch = "x86", target_arch = "x86_64")))]
pub fn apply_repulsion_2d_simd_f32_parallel(layout: &mut Layout) {
	#[cfg(target_arch = "x86")]
	use std::arch::x86::*;
	#[cfg(target_arch = "x86_64")]
	use std::arch::x86_64::*;

	let chunk_size = layout.settings.chunk_size.unwrap();
	let kr = layout.settings.kr;
	for chunk_iter in layout.iter_par_simd_nodes::<4>(chunk_size) {
		chunk_iter.for_each(|n1_iter| {
			for mut n1 in n1_iter {
				let n1_mass = n1.mass + 1.0f32;
				let n1_pos = unsafe {
					_mm256_set_ps(
						*n1.pos.get_unchecked(1),
						*n1.pos.get_unchecked(0),
						*n1.pos.get_unchecked(1),
						*n1.pos.get_unchecked(0),
						*n1.pos.get_unchecked(1),
						*n1.pos.get_unchecked(0),
						*n1.pos.get_unchecked(1),
						*n1.pos.get_unchecked(0),
					)
				};

				// This loop iterates on nodes by 4
				for n2 in &mut n1.n2_iter {
					unsafe {
						// [n2_x, n2_y, n3_x, n3_y, n4_x, n4_y, n5_x, n5_y]
						let n2345_pos = _mm256_loadu_ps(n2.pos);

						// [dx(n1,n2), dy(n1,n2), dx(n1,n3), dy(n1,n3), dx(n1,n4), dy(n1,n4), dx(n1,n5), dy(n1,n5)]
						let dxy = _mm256_sub_ps(n2345_pos, n1_pos);

						// ([dx(n1,n2)^2, dx(n1,n3)^2, dx(n1,n4)^2, dx(n1,n5)^2], [dy(n1,n2)^2, dy(n1,n3)^2, dy(n1,n4)^2, dy(n1,n5)^2])
						let (dx2, dy2): (__m128, __m128) =
							std::mem::transmute(_mm256_permutevar8x32_ps(
								_mm256_mul_ps(dxy, dxy),
								_mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0),
							));

						// [d(n1,n2), d(n1,n3), d(n1,n4), d(n1,n5)]
						let d2 = _mm_add_ps(dx2, dy2);
						// TODO maybe check zero

						let degs2 = n1_mass * (*n2.mass + 1.);
						let degs3 = n1_mass * (*n2.mass.add(1) + 1.);
						let degs4 = n1_mass * (*n2.mass.add(2) + 1.);
						let degs5 = n1_mass * (*n2.mass.add(3) + 1.);
						let f = _mm_mul_ps(
							_mm_div_ps(_mm_set_ps(degs5, degs4, degs3, degs2), d2),
							_mm_set_ps(kr, kr, kr, kr),
						);

						let (n1_speed_v, n2345_speed_v): (__m128, __m256) =
							(_mm_loadu_ps(n1.speed), _mm256_loadu_ps(n2.speed));

						let fxy = _mm256_mul_ps(
							dxy,
							_mm256_permutevar8x32_ps(
								_mm256_set_m128(f, f),
								_mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0),
							),
						);

						_mm256_storeu_ps(n2.speed, _mm256_add_ps(n2345_speed_v, fxy));
						_mm_storeu_ps(
							n1.speed,
							_mm_sub_ps(
								n1_speed_v,
								_mm_add_ps(
									_mm256_extractf128_ps(fxy, 1),
									_mm256_extractf128_ps(fxy, 0),
								),
							),
						);
					}
				}

				// Remaining iterations
				let layout = unsafe { n1.n2_iter.layout.0.as_mut() };
				for n2 in n1.n2_iter.ind
					..(n1.n2_iter.ind + (n1.n2_iter.ind - n1.ind + 1) / 4 * 4 + 4)
						.min(layout.masses.len())
				{
					let n2_pos = unsafe { layout.points.get_unchecked(n2) };

					let dx = unsafe { *n2_pos.get_unchecked(0) - *n1.pos.get_unchecked(0) };
					let dy = unsafe { *n2_pos.get_unchecked(1) - *n1.pos.get_unchecked(1) };

					let d2 = dx * dx + dy * dy;
					if d2.is_zero() {
						continue;
					}

					let f = (n1_mass * (unsafe { layout.masses.get_unchecked(n2) } + 1.)) / d2 * kr;

					let n2_speed = layout.speeds.get_mut(n2);
					let vx = f * dx;
					let vy = f * dy;
					unsafe { *n1.speed }.sub_assign(vx); // n1_speed[0] -= f * dx
					unsafe { *n1.speed.add(1) }.sub_assign(vy); // n1_speed[1] -= f * dy
					unsafe { n2_speed.get_unchecked_mut(0) }.add_assign(vx); // n2_speed[0] += f * dx
					unsafe { n2_speed.get_unchecked_mut(1) }.add_assign(vy); // n2_speed[1] += f * dy
				}
			}
		});
	}
}

pub fn apply_repulsion_3d(layout: &mut Layout) {
	for (n1, (n1_mass, n1_pos)) in layout.masses.iter().zip(layout.points.iter()).enumerate() {
		let mut n2_iter = layout.points.iter();
		let n1_mass = *n1_mass + Coord::one();
		for (n2, n2_pos) in (0..n1).zip(&mut n2_iter) {
			let dx = unsafe { *n2_pos.get_unchecked(0) - *n1_pos.get_unchecked(0) };
			let dy = unsafe { *n2_pos.get_unchecked(1) - *n1_pos.get_unchecked(1) };
			let dz = unsafe { *n2_pos.get_unchecked(2) - *n1_pos.get_unchecked(2) };

			let d2 = dx * dx + dy * dy + dz * dz;
			if d2.is_zero() {
				continue;
			}

			let f = n1_mass * (*unsafe { layout.masses.get_unchecked(n2) } + Coord::one()) / d2
				* layout.settings.kr;

			let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);
			unsafe { n1_speed.get_unchecked_mut(0) }.sub_assign(f * dx); // n1_speed[0] += f * dx
			unsafe { n1_speed.get_unchecked_mut(1) }.sub_assign(f * dy); // n1_speed[1] += f * dy
			unsafe { n1_speed.get_unchecked_mut(2) }.sub_assign(f * dz); // n1_speed[2] += f * dz
			unsafe { n2_speed.get_unchecked_mut(0) }.add_assign(f * dx); // n2_speed[0] -= f * dx
			unsafe { n2_speed.get_unchecked_mut(1) }.add_assign(f * dy); // n2_speed[1] -= f * dy
			unsafe { n2_speed.get_unchecked_mut(2) }.add_assign(f * dz); // n2_speed[2] -= f * dz
		}
	}
}

pub fn apply_repulsion_po(layout: &mut Layout) {
	let mut di = valloc(layout.settings.dimensions);
	let (node_size, krprime) = {
		layout
			.settings
			.prevent_overlapping
			.as_ref()
            // TODO - consider further `hint::unreachable_unchecked()` hinting as an optimization
			.unwrap()
	};
	for (n1, (n1_mass, n1_pos)) in layout.masses.iter().zip(layout.points.iter()).enumerate() {
		let mut n2_iter = layout.points.iter();
		let n1_mass = n1_mass.clone() + Coord::one();
		n2_iter.offset = (n1 + 1) * layout.settings.dimensions;
		for (n2, n2_pos) in (0..n1).zip(&mut n2_iter) {
			di.clone_from_slice(n2_pos);

			let d2 = di
				.iter_mut()
				.zip(n1_pos.iter())
				.map(|(di, n1_pos)| {
					*di -= n1_pos.clone();
					(*di).powi(2)
				})
				.sum::<Coord>();
			if d2.is_zero() {
				continue;
			}

			let d = d2.clone().sqrt();
			let dprime = d.clone() - node_size.clone();

			let f = n1_mass.clone()
				* (unsafe { layout.masses.get_unchecked(n2) }.clone() + Coord::one())
				/ d2 * if dprime.is_sign_positive() {
				layout.settings.kr.clone() / dprime
			} else {
				krprime.clone()
			};

			let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);
			izip!(n1_speed.iter_mut(), n2_speed.iter_mut(), di.iter()).for_each(
				|(n1_speed, n2_speed, di)| {
					let s = f.clone() * di.clone();
					*n1_speed -= s.clone();
					*n2_speed += s;
				},
			);
		}
	}
}

#[cfg(feature = "barnes_hut")]
pub fn apply_repulsion_bh_2d(layout: &mut Layout<f64>) {
	let particles: Vec<nbody_barnes_hut::particle_2d::Particle2D> = layout
		.points
		.iter()
		.zip(layout.masses.iter())
		.map(|(point, mass)| nbody_barnes_hut::particle_2d::Particle2D {
			position: nbody_barnes_hut::vector_2d::Vector2D {
				x: point[0],
				y: point[1],
			},
			mass: mass + 1.,
		})
		.collect();
	let tree = nbody_barnes_hut::barnes_hut_2d::QuadTree::new(
		&particles
			.iter()
			.collect::<Vec<&nbody_barnes_hut::particle_2d::Particle2D>>(),
		layout.settings.barnes_hut.unwrap(),
	);
	let kr = layout.settings.kr;

	izip!(
		particles.into_iter(),
		layout.speeds.iter_mut(),
		layout.masses.iter()
	)
	.for_each(|(particle, speed, mass)| {
		let nbody_barnes_hut::vector_2d::Vector2D { x, y } =
			tree.calc_forces_on_particle(particle.position, mass + 1., |d2, m1, dv, m2| {
				unsafe {
					std::intrinsics::fmul_fast(
						std::intrinsics::fmul_fast(m2, std::intrinsics::fdiv_fast(m1, d2.sqrt())),
						kr,
					) * dv
				}
				//m2 as f64 * m1 / d2.sqrt() * kr * dv
			});
		speed[0] -= x;
		speed[1] -= y;
	});
}

#[cfg(feature = "barnes_hut")]
pub fn apply_repulsion_bh_2d_po(layout: &mut Layout<f64>) {
	let particles: Vec<nbody_barnes_hut::particle_2d::Particle2D> = layout
		.points
		.iter()
		.zip(layout.masses.iter())
		.map(|(point, mass)| nbody_barnes_hut::particle_2d::Particle2D {
			position: nbody_barnes_hut::vector_2d::Vector2D {
				x: point[0],
				y: point[1],
			},
			mass: mass + 1.,
		})
		.collect();
	let tree = nbody_barnes_hut::barnes_hut_2d::QuadTree::new(
		&particles
			.iter()
			.collect::<Vec<&nbody_barnes_hut::particle_2d::Particle2D>>(),
		layout.settings.barnes_hut.unwrap(),
	);
	let kr = layout.settings.kr;
	let (node_size, krprime) = unsafe { layout.settings.prevent_overlapping.unwrap_unchecked() };
	izip!(
		particles.into_iter(),
		layout.speeds.iter_mut(),
		layout.masses.iter()
	)
	.for_each(|(particle, speed, mass)| {
		let nbody_barnes_hut::vector_2d::Vector2D { x, y } =
			tree.calc_forces_on_particle(particle.position, mass + 1., |d2, m1, dv, m2| {
				let d = d2.sqrt();
				let dprime = d - node_size;
				(if dprime.positive() {
					kr / dprime
				} else if dprime.is_zero() {
					return nbody_barnes_hut::vector_2d::Vector2D { x: 0.0, y: 0.0 };
				} else {
					krprime
				}) * m1 * m2 / d * dv
			});
		speed[0] -= x;
		speed[1] -= y;
	});
}

#[cfg(feature = "barnes_hut")]
pub fn apply_repulsion_bh_3d(layout: &mut Layout<f64>) {
	let particles: Vec<nbody_barnes_hut::particle_3d::Particle3D> = layout
		.points
		.iter()
		.zip(layout.masses.iter())
		.map(|(point, mass)| nbody_barnes_hut::particle_3d::Particle3D {
			position: nbody_barnes_hut::vector_3d::Vector3D {
				x: point[0],
				y: point[1],
				z: point[2],
			},
			mass: mass + 1.,
		})
		.collect();
	let tree = nbody_barnes_hut::barnes_hut_3d::OctTree::new(
		&particles
			.iter()
			.collect::<Vec<&nbody_barnes_hut::particle_3d::Particle3D>>(),
		layout.settings.barnes_hut.unwrap(),
	);
	let kr = layout.settings.kr;
	izip!(
		particles.into_iter(),
		layout.speeds.iter_mut(),
		layout.masses.iter()
	)
	.for_each(|(particle, speed, mass)| {
		let nbody_barnes_hut::vector_3d::Vector3D { x, y, z } =
			tree.calc_forces_on_particle(particle.position, mass + 1., |d2, m1, dv, m2| {
				m2 * m1 / d2.sqrt() * kr * dv
			});
		speed[0] -= x;
		speed[1] -= y;
		speed[2] -= z;
	});
}

#[cfg(feature = "barnes_hut")]
pub fn apply_repulsion_bh_3d_po(layout: &mut Layout<f64>) {
	let particles: Vec<nbody_barnes_hut::particle_3d::Particle3D> = layout
		.points
		.iter()
		.zip(layout.masses.iter())
		.map(|(point, mass)| nbody_barnes_hut::particle_3d::Particle3D {
			position: nbody_barnes_hut::vector_3d::Vector3D {
				x: point[0],
				y: point[1],
				z: point[2],
			},
			mass: mass + 1.,
		})
		.collect();
	let tree = nbody_barnes_hut::barnes_hut_3d::OctTree::new(
		&particles
			.iter()
			.collect::<Vec<&nbody_barnes_hut::particle_3d::Particle3D>>(),
		layout.settings.barnes_hut.unwrap(),
	);
	let kr = layout.settings.kr;
	let (node_size, krprime) = unsafe { layout.settings.prevent_overlapping.unwrap_unchecked() };
	izip!(
		particles.into_iter(),
		layout.speeds.iter_mut(),
		layout.masses.iter()
	)
	.for_each(|(particle, speed, mass)| {
		let nbody_barnes_hut::vector_3d::Vector3D { x, y, z } =
			tree.calc_forces_on_particle(particle.position, mass + 1., |d2, m1, dv, m2| {
				let d = d2.sqrt();
				let dprime = d - node_size;
				(if dprime.positive() {
					kr / dprime
				} else if dprime.is_zero() {
					return nbody_barnes_hut::vector_3d::Vector3D {
						x: 0.0,
						y: 0.0,
						z: 0.0,
					};
				} else {
					krprime
				}) * m1 * m2 / d * dv
			});
		speed[0] -= x;
		speed[1] -= y;
		speed[2] -= z;
	});
}
