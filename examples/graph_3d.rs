#![allow(dead_code)]

use forceatlas2::*;
use rand::Rng;
use std::{f32::consts::PI, fs::File, io::Read};
use stl::*;

const EDGES: usize = 5000;
const NODES: usize = 1000;

const ITERATIONS: u32 = 10;

const NODE_RADIUS: f32 = 10.0;
const EDGE_RADIUS: f32 = 1.0;
const SCALE: f32 = 0.001;

#[inline]
fn tri(normal: [f32; 3], v1: [f32; 3], v2: [f32; 3], v3: [f32; 3]) -> Triangle {
	Triangle {
		normal,
		v1,
		v2,
		v3,
		attr_byte_count: 0,
	}
}

fn octahedron(r: f32) -> Vec<Triangle> {
	vec![
		tri([1., 1., 1.], [0., 0., r], [0., r, 0.], [r, 0., 0.]),
		tri([1., 1., -1.], [0., 0., -r], [0., r, 0.], [r, 0., 0.]),
		tri([1., -1., 1.], [0., 0., r], [0., -r, 0.], [r, 0., 0.]),
		tri([1., -1., -1.], [0., 0., -r], [0., -r, 0.], [r, 0., 0.]),
		tri([-1., 1., 1.], [0., 0., r], [0., r, 0.], [-r, 0., 0.]),
		tri([-1., 1., -1.], [0., 0., -r], [0., r, 0.], [-r, 0., 0.]),
		tri([-1., -1., 1.], [0., 0., r], [0., -r, 0.], [-r, 0., 0.]),
		tri([-1., -1., -1.], [0., 0., -r], [0., -r, 0.], [-r, 0., 0.]),
	]
}

fn middle(v1: [f32; 3], v2: [f32; 3]) -> [f32; 3] {
	[
		(v1[0] + v2[0]) / 2.0,
		(v1[1] + v2[1]) / 2.0,
		(v1[2] + v2[2]) / 2.0,
	]
}

fn ray_normal(v1: [f32; 3], v2: [f32; 3], v3: [f32; 3]) -> [f32; 3] {
	[
		(v1[0] + v2[0] + v3[0]) / 3.0,
		(v1[1] + v2[1] + v3[1]) / 3.0,
		(v1[2] + v2[2] + v3[2]) / 3.0,
	]
}

fn project(v: [f32; 3], r: f32) -> [f32; 3] {
	let f = r / (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt();
	[v[0] * f, v[1] * f, v[2] * f]
}

fn ico_sphere(r: f32, iterations: usize) -> Vec<Triangle> {
	let mut sphere = octahedron(r);
	let mut sphere2 = Vec::<Triangle>::new();

	for _ in 0..iterations {
		for triangle in sphere.drain(..) {
			let v1 = project(middle(triangle.v3, triangle.v2), r);
			let v2 = project(middle(triangle.v3, triangle.v1), r);
			let v3 = project(middle(triangle.v2, triangle.v1), r);

			sphere2.push(tri(ray_normal(v1, v2, v3), v1, v2, v3));
			sphere2.push(tri(ray_normal(v1, v2, triangle.v3), v1, v2, triangle.v3));
			sphere2.push(tri(ray_normal(v1, triangle.v2, v3), v1, triangle.v2, v3));
			sphere2.push(tri(ray_normal(triangle.v1, v2, v3), triangle.v1, v2, v3));
		}
		sphere.append(&mut sphere2);
	}

	sphere
}

fn translate(solid: &mut Vec<Triangle>, offset: [f32; 3]) {
	for triangle in solid {
		triangle.v1[0] += offset[0];
		triangle.v1[1] += offset[1];
		triangle.v1[2] += offset[2];
		triangle.v2[0] += offset[0];
		triangle.v2[1] += offset[1];
		triangle.v2[2] += offset[2];
		triangle.v3[0] += offset[0];
		triangle.v3[1] += offset[1];
		triangle.v3[2] += offset[2];
	}
}

/// Allocate Vec without initializing
#[inline]
fn valloc<T>(n: usize) -> Vec<T> {
	let mut v = Vec::with_capacity(n);
	unsafe {
		v.set_len(n);
	}
	v
}

fn clone_solid(solid: &Vec<Triangle>) -> Vec<Triangle> {
	let mut new_solid = valloc(solid.len());
	for i in 0..solid.len() {
		let triangle = unsafe { solid.get_unchecked(i) };
		new_solid[i] = Triangle {
			normal: triangle.normal,
			v1: triangle.v1,
			v2: triangle.v2,
			v3: triangle.v3,
			attr_byte_count: triangle.attr_byte_count,
		};
	}
	new_solid
}

fn rotate_point(point: &mut [f32; 3], angle: [f32; 3]) {
	point[0] = point[0] * angle[2].cos() - point[1] * angle[2].sin();
	point[1] = point[0] * angle[2].sin() - point[1] * angle[2].cos();
	point[2] = point[2] * angle[0].cos() - point[1] * angle[0].sin();
	point[1] = point[2] * angle[0].sin() - point[1] * angle[0].cos();
	point[0] = point[0] * angle[1].cos() - point[2] * angle[1].sin();
	point[2] = point[0] * angle[1].sin() - point[2] * angle[1].cos();
}

fn rotate(solid: &mut Vec<Triangle>, angle: [f32; 3]) {
	for triangle in solid {
		rotate_point(&mut triangle.normal, angle);
		rotate_point(&mut triangle.v1, angle);
		rotate_point(&mut triangle.v2, angle);
		rotate_point(&mut triangle.v3, angle);
	}
}

fn cylinder(r: f32, h: f32, n: usize) -> Vec<Triangle> {
	let mut cylinder = valloc(4 * n);

	let ai = 2.0 * PI / n as f32;
	let mut lv = [r, 0.0, 0.0];
	let mut j = 0;
	for i in 1..=n {
		let a = i as f32 * ai;
		let v = [r * a.cos(), r * a.sin(), 0.0];
		cylinder[j] = tri([0.0, 0.0, -1.0], [0.0, 0.0, 0.0], lv, v);
		cylinder[j + 1] = tri(
			[0.0, 0.0, 1.0],
			[0.0, 0.0, h],
			[lv[0], lv[1], h],
			[v[0], v[1], h],
		);
		let normal = middle(lv, v);
		cylinder[j + 2] = tri(normal, lv, v, [v[0], v[1], h]);
		cylinder[j + 3] = tri(normal, lv, [v[0], v[1], h], [lv[0], lv[1], h]);
		lv = v;
		j += 4
	}

	cylinder
}

fn scale_point(point: &mut [f32; 3], factor: [f32; 3]) {
	point[0] *= factor[0];
	point[1] *= factor[1];
	point[2] *= factor[2];
}

fn scale(solid: &mut Vec<Triangle>, factor: [f32; 3]) {
	for triangle in solid {
		scale_point(&mut triangle.v1, factor);
		scale_point(&mut triangle.v2, factor);
		scale_point(&mut triangle.v3, factor);
	}
}

fn uniscale(solid: &mut Vec<Triangle>, factor: f32) {
	for triangle in solid {
		triangle.v1[0] *= factor;
		triangle.v1[1] *= factor;
		triangle.v1[2] *= factor;
		triangle.v2[0] *= factor;
		triangle.v2[1] *= factor;
		triangle.v2[2] *= factor;
		triangle.v3[0] *= factor;
		triangle.v3[1] *= factor;
		triangle.v3[2] *= factor;
	}
}

fn distance(v1: [f32; 3], v2: [f32; 3]) -> f32 {
	((v2[0] - v1[0]).powi(2) + (v2[1] - v1[1]).powi(2) + (v2[2] - v1[2]).powi(2)).sqrt()
}

fn main() {
	// -------- Generate graph
	let (edges, nb_nodes) = if let Some(path) = std::env::args().skip(1).next() {
		let mut file = std::fs::File::open(path).expect("Cannot open file");
		let mut contents = String::new();
		file.read_to_string(&mut contents).unwrap();
		let mut max_node = 0usize;
		(
			contents
				.lines()
				.skip(1)
				.map(|line| {
					let mut columns = line.split(';');
					let edge = (
						columns.next().unwrap().parse().unwrap(),
						columns.next().unwrap().parse().unwrap(),
					);
					if edge.0 > max_node {
						max_node = edge.0;
					}
					if edge.1 > max_node {
						max_node = edge.1;
					}
					edge
				})
				.collect::<Vec<(usize, usize)>>(),
			max_node + 1,
		)
	} else {
		let mut rng = rand::thread_rng();
		(
			(0..EDGES)
				.map(|_| {
					let edge = (rng.gen_range(0usize..NODES), rng.gen_range(0..NODES));
					if edge.0 == edge.1 {
						(edge.0, (edge.1 + 1) % NODES)
					} else {
						edge
					}
				})
				.collect(),
			NODES,
		)
	};

	// Compute degrees
	let mut degrees = vec![0usize; nb_nodes];
	for (n1, n2) in edges.iter() {
		degrees[*n1] += 1;
		degrees[*n2] += 1;
	}

	// Compute layout
	let mut layout = Layout::<f64>::from_graph(
		edges,
		nb_nodes,
		Settings {
			dimensions: 3,
			dissuade_hubs: true,
			ka: 0.5,
			kg: 10000.0, //1.0
			kr: 100.0,   //0.1
			lin_log: false,
			prevent_overlapping: Some((NODE_RADIUS as f64, 100.0)),
			strong_gravity: false,
			barnes_hut: Some(0.5),
		},
	);

	for _ in 0..ITERATIONS {
		layout.iteration();
	}

	// -------- Generate Solid
	let sphere = ico_sphere(NODE_RADIUS, 2);
	//let cylinder = cylinder(EDGE_RADIUS, 1.0, 6);
	let mut solid = Vec::with_capacity(
		sphere.len() * nb_nodes, /* + cylinder.len() * edges.len()*/
	);

	for (i, node) in layout.points.iter().enumerate() {
		let mut new_sphere = clone_solid(&sphere);
		uniscale(&mut new_sphere, degrees[i] as f32);
		translate(
			&mut new_sphere,
			dbg!([
				node[0] as f32 * SCALE,
				node[1] as f32 * SCALE,
				node[2] as f32 * SCALE
			]),
		);
		solid.append(&mut new_sphere);
	}

	/*for (n1, n2) in layout.edges {
		let node1 = layout.points.get(n1);
		let node2 = layout.points.get(n2);
		let node1 = [node1[0] as f32, node1[1] as f32, node1[2] as f32];
		let node2 = [node2[0] as f32, node2[1] as f32, node2[2] as f32];

		let mut new_cylinder = clone_solid(&cylinder);
		scale(&mut new_cylinder, [1.0, 1.0, distance(node1, node2)]);
		rotate(
			&mut new_cylinder,
			[
				PI - (node2[1] - node1[1]).atan2(node2[2] - node1[2]),
				PI/2.0 - (node2[2] - node1[2]).atan2(node2[0] - node1[0]),
				0.0,//(node2[1] - node1[1]).atan2(node2[0] - node1[0]),
			],
		);
		translate(
			&mut new_cylinder,
			[
				node1[0],
				node1[1],
				node1[2],
			],
		);
		dbg!(new_cylinder.get(0).unwrap().v1);
		dbg!(new_cylinder.get(0).unwrap().v2);
		dbg!(new_cylinder.get(0).unwrap().v3);
		solid.append(&mut new_cylinder);
	}*/

	let mut out = File::create("target/graph.stl").unwrap();
	write_stl(
		&mut out,
		&BinaryStlFile {
			header: BinaryStlHeader {
				header: [0; 80],
				num_triangles: solid.len() as u32,
			},
			triangles: solid,
		},
	)
	.unwrap();
}
