use forceatlas2::*;
use plotters::prelude::*;
use rand::Rng;

const EDGES: usize = 50_000;
const NODES: usize = 5_000;

const SIZE: (u32, u32) = (1024, 1024);

const ITERATIONS: u32 = 50;
const ANIM_MODE: bool = false;

fn main() {
	let mut rng = rand::thread_rng();
	eprintln!("Generating graph...");
	let edges = (0..EDGES)
		.map(|_| {
			let edge = (rng.gen_range(0usize..NODES), rng.gen_range(0..NODES));
			if edge.0 == edge.1 {
				(edge.0, (edge.1 + 1) % NODES)
			} else {
				edge
			}
		})
		.collect();
	let mut layout = Layout::<f64>::from_graph(
		edges,
		NODES,
		Settings {
			dimensions: 2,
			dissuade_hubs: false,
			ka: 0.5,
			kg: 1.0,
			kr: 0.1,
			lin_log: false,
			prevent_overlapping: None,
			strong_gravity: false,
			barnes_hut: None,
		},
	);

	eprintln!("Computing layout...");
	for i in 0..ITERATIONS {
		if ANIM_MODE {
			draw_graph(&layout, i);
		}
		println!("{}/{}", i, ITERATIONS);
		layout.iteration();
	}
	draw_graph(&layout, ITERATIONS);
}

fn draw_graph(layout: &Layout<f64>, iteration: u32) {
	let mut min_v = layout.points.get_clone(0);
	let mut max_v = min_v.clone();
	let min = min_v.as_mut_slice();
	let max = max_v.as_mut_slice();
	for pos in layout.points.iter() {
		if pos[0] < min[0] {
			min[0] = pos[0];
		}
		if pos[1] < min[1] {
			min[1] = pos[1];
		}
		if pos[0] > max[0] {
			max[0] = pos[0];
		}
		if pos[1] > max[1] {
			max[1] = pos[1];
		}
	}
	let graph_size = (max[0] - min[0], max[1] - min[1]);
	let factor = {
		let factors = (
			f64::from(SIZE.0) / graph_size.0,
			f64::from(SIZE.1) / graph_size.1,
		);
		if factors.0 > factors.1 {
			min[0] -= (f64::from(SIZE.0) / factors.1 - graph_size.0) / 2.0;
			factors.1
		} else {
			min[1] -= (f64::from(SIZE.1) / factors.0 - graph_size.1) / 2.0;
			factors.0
		}
	};
	println!("size:  {:?}", graph_size);
	println!("scale: {}", factor);

	let path = if ANIM_MODE {
		format!("target/graph-{}.png", iteration)
	} else {
		"target/graph.png".into()
	};
	let root = BitMapBackend::new(&path, SIZE).into_drawing_area();
	root.fill(&WHITE).unwrap();

	if !ANIM_MODE {
		for (h1, h2) in layout.edges.iter() {
			root.draw(&PathElement::new(
				vec![
					{
						let pos = layout.points.get(*h1);
						unsafe {
							(
								((pos[0] - min[0]) * factor).to_int_unchecked::<i32>(),
								((pos[1] - min[1]) * factor).to_int_unchecked::<i32>(),
							)
						}
					},
					{
						let pos = layout.points.get(*h2);
						unsafe {
							(
								((pos[0] - min[0]) * factor).to_int_unchecked::<i32>(),
								((pos[1] - min[1]) * factor).to_int_unchecked::<i32>(),
							)
						}
					},
				],
				Into::<ShapeStyle>::into(&BLACK).filled(),
			))
			.unwrap();
		}
	}

	for pos in layout.points.iter() {
		root.draw(&Circle::new(
			unsafe {
				(
					((pos[0] - min[0]) * factor).to_int_unchecked::<i32>(),
					((pos[1] - min[1]) * factor).to_int_unchecked::<i32>(),
				)
			},
			2,
			Into::<ShapeStyle>::into(&RED).filled(),
		))
		.unwrap();
	}
}
