use crate::{layout::Layout, util::*};

use itertools::izip;

pub fn apply_gravity<T: Coord + std::fmt::Debug>(layout: &mut Layout<T>) {
	for (n, pos, speed) in izip!(
		layout.nodes.iter(),
		layout.points.iter(),
		layout.speeds.iter_mut()
	) {
		let d = norm(pos);
		if d.is_zero() {
			continue;
		}
		let f = T::from(n.degree + 1) * layout.settings.kg.clone() / d;
		for (speed, pos) in speed.iter_mut().zip(pos.iter()) {
			*speed -= f.clone() * pos.clone();
		}
	}
}

pub fn apply_gravity_sg<T: Coord + std::fmt::Debug>(layout: &mut Layout<T>) {
	for (n, pos, speed) in izip!(
		layout.nodes.iter(),
		layout.points.iter(),
		layout.speeds.iter_mut()
	) {
		let f = T::from(n.degree + 1) * layout.settings.kg.clone();
		for (speed, pos) in speed.iter_mut().zip(pos.iter()) {
			*speed -= f.clone() * pos.clone();
		}
	}
}
