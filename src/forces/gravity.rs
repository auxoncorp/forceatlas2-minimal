use crate::{layout::Layout, util::*};

use itertools::izip;
use num_traits::{One, Zero};

pub fn apply_gravity(layout: &mut Layout) {
	for (mass, pos, speed) in izip!(
		layout.masses.iter(),
		layout.points.iter(),
		layout.speeds.iter_mut()
	) {
		let d = norm(pos);
		if d.is_zero() {
			continue;
		}
		let f = (mass.clone() + Coord::one()) * layout.settings.kg.clone() / d;
		for (speed, pos) in speed.iter_mut().zip(pos.iter()) {
			*speed -= f.clone() * pos.clone();
		}
	}
}

pub fn apply_gravity_sg(layout: &mut Layout) {
	for (mass, pos, speed) in izip!(
		layout.masses.iter(),
		layout.points.iter(),
		layout.speeds.iter_mut()
	) {
		let f = (mass.clone() + Coord::one()) * layout.settings.kg.clone();
		for (speed, pos) in speed.iter_mut().zip(pos.iter()) {
			*speed -= f.clone() * pos.clone();
		}
	}
}
