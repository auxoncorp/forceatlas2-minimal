use crate::{layout::Layout, util::*};

use itertools::izip;
use num_traits::{Zero};
use std::ops::{AddAssign, SubAssign};

pub fn apply_attraction(layout: &mut Layout) {
	let mut di_v = valloc(layout.settings.dimensions);
	let di = di_v.as_mut_slice();
	for (n1, n2) in layout.edges.iter() {
		let (n1, n2) = (*n1, *n2);
		let n1_pos = layout.points.get(n1);
		layout.points.get_clone_slice(n2, di);

		let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);
		for i in 0..layout.settings.dimensions {
			let di = unsafe { di.get_unchecked_mut(i) };
			let n1_speed = unsafe { n1_speed.get_unchecked_mut(i) };
			let n2_speed = unsafe { n2_speed.get_unchecked_mut(i) };
			let n1_pos = unsafe { n1_pos.get_unchecked(i) };

			*di -= n1_pos.clone();
			*di *= layout.settings.ka.clone();
			*n1_speed += di.clone();
			*n2_speed -= di.clone();
		}
	}
}

pub fn apply_attraction_2d(layout: &mut Layout) {
	for (n1, n2) in layout.edges.iter() {
		let (n1, n2) = (*n1, *n2);

		let n1_pos = layout.points.get(n1);
		let n2_pos = layout.points.get(n2);

		let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);

		let dx =
			unsafe { *n2_pos.get_unchecked(0) - *n1_pos.get_unchecked(0) } * layout.settings.ka;
		let dy =
			unsafe { *n2_pos.get_unchecked(1) - *n1_pos.get_unchecked(1) } * layout.settings.ka;

		unsafe { n1_speed.get_unchecked_mut(0) }.add_assign(dx);
		unsafe { n1_speed.get_unchecked_mut(1) }.add_assign(dy);
		unsafe { n2_speed.get_unchecked_mut(0) }.sub_assign(dx);
		unsafe { n2_speed.get_unchecked_mut(1) }.sub_assign(dy);
	}
}

pub fn apply_attraction_3d(layout: &mut Layout) {
	for (n1, n2) in layout.edges.iter() {
		let (n1, n2) = (*n1, *n2);

		let n1_pos = layout.points.get(n1);
		let n2_pos = layout.points.get(n2);

		let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);

		let dx =
			unsafe { *n2_pos.get_unchecked(0) - *n1_pos.get_unchecked(0) } * layout.settings.ka;
		let dy =
			unsafe { *n2_pos.get_unchecked(1) - *n1_pos.get_unchecked(1) } * layout.settings.ka;
		let dz =
			unsafe { *n2_pos.get_unchecked(2) - *n1_pos.get_unchecked(2) } * layout.settings.ka;

		unsafe { n1_speed.get_unchecked_mut(0) }.add_assign(dx);
		unsafe { n1_speed.get_unchecked_mut(1) }.add_assign(dy);
		unsafe { n1_speed.get_unchecked_mut(2) }.add_assign(dz);
		unsafe { n2_speed.get_unchecked_mut(0) }.sub_assign(dx);
		unsafe { n2_speed.get_unchecked_mut(1) }.sub_assign(dy);
		unsafe { n2_speed.get_unchecked_mut(2) }.sub_assign(dz);
	}
}

pub fn apply_attraction_dh(layout: &mut Layout) {
	for (n1, n2) in layout.edges.iter() {
		let n1_speed = layout.speeds.get_mut(*n1);
		let n1_pos = layout.points.get(*n1);
		let mut di_v = layout.points.get_clone(*n2);
		let di = di_v.as_mut_slice();
		let n1_mass = layout.masses.get(*n1).unwrap().clone();
		for (n1_speed, n1_pos, di) in izip!(n1_speed, n1_pos, di.iter_mut()) {
			*di -= n1_pos.clone();
			*di /= n1_mass.clone();
			*di *= layout.settings.ka.clone();
			*n1_speed += di.clone();
		}
		let n2_speed = layout.speeds.get_mut(*n2);
		for i in 0usize..layout.settings.dimensions {
			n2_speed[i] -= di[i].clone();
		}
	}
}

pub fn apply_attraction_log(layout: &mut Layout) {
	for (n1, n2) in layout.edges.iter() {
		let mut d = Coord::zero();
		let mut di_v = layout.points.get_clone(*n2);
		let di = di_v.as_mut_slice();
		for (di, n1) in di.iter_mut().zip(layout.points.get(*n1)) {
			*di -= n1.clone();
			d += (*di).powi(2);
		}
		if d.is_zero() {
			continue;
		}
		d = d.sqrt();

		let f = d.clone().ln_1p() / d * layout.settings.ka.clone();

		let n1_speed = layout.speeds.get_mut(*n1);
		for i in 0usize..layout.settings.dimensions {
			n1_speed[i] += f.clone() * di[i].clone();
		}
		let n2_speed = layout.speeds.get_mut(*n2);
		for i in 0usize..layout.settings.dimensions {
			n2_speed[i] -= f.clone() * di[i].clone();
		}
	}
}

pub fn apply_attraction_dh_log(layout: &mut Layout) {
	for (n1, n2) in layout.edges.iter() {
		let mut d = Coord::zero();
		let mut di_v = layout.points.get_clone(*n2);
		let di = di_v.as_mut_slice();
		for (di, n1) in di.iter_mut().zip(layout.points.get(*n1)) {
			*di -= n1.clone();
			d += (*di).powi(2);
		}
		if d.is_zero() {
			continue;
		}
		d = d.sqrt();

		let n1_mass = layout.masses.get(*n1).unwrap().clone();
		let f = d.clone().ln_1p() / d / n1_mass * layout.settings.ka.clone();

		let n1_speed = layout.speeds.get_mut(*n1);
		for i in 0usize..layout.settings.dimensions {
			n1_speed[i] += f.clone() * di[i].clone();
		}
		let n2_speed = layout.speeds.get_mut(*n2);
		for i in 0usize..layout.settings.dimensions {
			n2_speed[i] -= f.clone() * di[i].clone();
		}
	}
}

pub fn apply_attraction_po(layout: &mut Layout) {
	let node_size = &layout.settings.prevent_overlapping.as_ref().unwrap().0;
	for (n1, n2) in layout.edges.iter() {
		let mut d = Coord::zero();
		let n1_pos = layout.points.get(*n1);
		let mut di_v = layout.points.get_clone(*n2);
		let di = di_v.as_mut_slice();
		for i in 0usize..layout.settings.dimensions {
			di[i] -= n1_pos[i].clone();
			d += (di[i]).powi(2);
		}
		d = d.sqrt();

		let dprime = d.clone() - node_size.clone();
		if !dprime.is_sign_positive() {
			continue;
		}
		let f = dprime / d * layout.settings.ka.clone();

		let n1_speed = layout.speeds.get_mut(*n1);
		for i in 0usize..layout.settings.dimensions {
			n1_speed[i] += f.clone() * di[i].clone();
		}
		let n2_speed = layout.speeds.get_mut(*n2);
		for i in 0usize..layout.settings.dimensions {
			n2_speed[i] -= f.clone() * di[i].clone();
		}
	}
}

pub fn apply_attraction_dh_po(layout: &mut Layout) {
	let node_size = &layout.settings.prevent_overlapping.as_ref().unwrap().0;
	for (n1, n2) in layout.edges.iter() {
		let mut d = Coord::zero();
		let n1_pos = layout.points.get(*n1);
		let mut di_v = layout.points.get_clone(*n2);
		let di = di_v.as_mut_slice();
		for i in 0usize..layout.settings.dimensions {
			di[i] -= n1_pos[i].clone();
			d += (di[i]).powi(2);
		}
		d = d.sqrt();

		let dprime = d.clone() - node_size.clone();
		if !dprime.is_sign_positive() {
			dbg!(dprime);
			continue;
		}
		let n1_mass = layout.masses.get(*n1).unwrap().clone();
		let f = dprime / d / n1_mass * layout.settings.ka.clone();

		let n1_speed = layout.speeds.get_mut(*n1);
		for i in 0usize..layout.settings.dimensions {
			n1_speed[i] += f.clone() * di[i].clone();
		}
		let n2_speed = layout.speeds.get_mut(*n2);
		for i in 0usize..layout.settings.dimensions {
			n2_speed[i] -= f.clone() * di[i].clone();
		}
	}
}

pub fn apply_attraction_log_po(layout: &mut Layout) {
	let node_size = &layout.settings.prevent_overlapping.as_ref().unwrap().0;
	for (n1, n2) in layout.edges.iter() {
		let mut d = Coord::zero();
		let n1_pos = layout.points.get(*n1);
		let mut di_v = layout.points.get_clone(*n2);
		let di = di_v.as_mut_slice();
		for i in 0usize..layout.settings.dimensions {
			di[i] -= n1_pos[i].clone();
			d += (di[i]).powi(2);
		}
		d = d.sqrt();

		let dprime = d - node_size.clone();
		if !dprime.is_sign_positive() {
			continue;
		}
		let f = dprime.clone().ln_1p() / dprime * layout.settings.ka.clone();

		let n1_speed = layout.speeds.get_mut(*n1);
		for i in 0usize..layout.settings.dimensions {
			n1_speed[i] += f.clone() * di[i].clone();
		}
		let n2_speed = layout.speeds.get_mut(*n2);
		for i in 0usize..layout.settings.dimensions {
			n2_speed[i] -= f.clone() * di[i].clone();
		}
	}
}

pub fn apply_attraction_dh_log_po(layout: &mut Layout) {
	let node_size = &layout.settings.prevent_overlapping.as_ref().unwrap().0;
	for (n1, n2) in layout.edges.iter() {
		let mut d = Coord::zero();
		let n1_pos = layout.points.get(*n1);
		let mut di_v = layout.points.get_clone(*n2);
		let di = di_v.as_mut_slice();
		for i in 0usize..layout.settings.dimensions {
			di[i] -= n1_pos[i].clone();
			d += (di[i]).powi(2);
		}
		d = d.sqrt();

		let dprime = d - node_size.clone();
		if !dprime.is_sign_positive() {
			continue;
		}
		let n1_mass = layout.masses.get(*n1).unwrap().clone();
		let f = dprime.clone().ln_1p() / dprime / n1_mass * layout.settings.ka.clone();

		let n1_speed = layout.speeds.get_mut(*n1);
		for i in 0usize..layout.settings.dimensions {
			n1_speed[i] += f.clone() * di[i].clone();
		}
		let n2_speed = layout.speeds.get_mut(*n2);
		for i in 0usize..layout.settings.dimensions {
			n2_speed[i] -= f.clone() * di[i].clone();
		}
	}
}
