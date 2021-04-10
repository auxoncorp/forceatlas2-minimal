# ForceAtlas2 Rust

Very fast implementation of [ForceAtlas2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4051631/) &#8211; force-directed Continuous Graph Layout Algorithm for Handy Network Visualization (i.e. position the nodes of a n-dimension graph for drawing it more human-readably)

![Example graph spacialized with ForceAtlas2-rs](https://txmn.tk/img/wot-fa2rs.png)

## Optimization

The implementations used depend on the type and the parameters. The most optimized is a `Copy` type with `prevent_overlapping` and `barnes_hut` disabled, in 2D or 3D. Some specializations are not implemented yet. `x86` and `x86_64` processors with support to `avx2` use SIMD to compute faster on the `f64` or `f32`, 2D, `prevent_overlapping` disabled and `barnes_hut` disabled case.

TL;DR If you want best performance, use the following:
* CPU: `x86` or `x86_64` with `avx2`
* `Layout<f32>` (`Layout<f64>` should be ~2 times slower)
* `Settings::prevent_overlapping: None`
* `Settings::barnes_hut: false` (or just don't use this feature)
* `RUSTFLAGS='-C target-feature=+avx2'`

Use the `barnes_hut` feature to turn repulsion from O(n^2) to O(n×log(n)) (only for 2D/3D and `f64`/`f32`). However, some optimizations like SIMD are not available with Barnes-Hut.

## Examples

[Install Rustup](https://rustup.rs/) and switch to nightly:

    rustup toolchain install nightly && rustup default nightly

A packet may be needed to draw graph:

    sudo apt install libfreetype6-dev

Clone repository:

    git clone https://framagit.org/ZettaScript/forceatlas2-rs && cd forceatlas2-rs

Build example: (`examples/wot.csv` file lists the edges of a directed graph, in two columns)

    RUSTFLAGS='-C target-feature=+avx2' cargo build --release --example csv_import
    ./target/release/examples/csv_import examples/wot.csv

Output images are in `target` directory.

## Comparison

Python (forceatlas2, fa2) and JS (sigma.js) implementations are slow.

Java implementation (Gephi) does not use SIMD.

Julia implementation (Anim-Wotmap) beats them all and uses SIMD, but is still slower than this one.

If you know any faster comparable force-directed layout implementation, please let me know.

## Bindings

There is a binding for use in Python, [fa2rs](https://framagit.org/ZettaScript/fa2rs-py).

## License

GNU AGPL v3, CopyLeft 2020-2021 Pascal Engélibert

Implementation details inspired by:
* [python-forceatlas2](https://code.launchpad.net/forceatlas2-python) (GNU GPL v3, CopyLeft 2016 Max Shinn)
* [python-fa2](https://github.com/bhargavchippada/forceatlas2) (GNU GPL v3, CopyLeft 2017 Bhargav Chippada)
* [Gephi](https://github.com/gephi/gephi/tree/master/modules/LayoutPlugin/src/main/java/org/gephi/layout/plugin/forceAtlas2) (GNU GPL v3 / CDDL 1.0, CopyLeft 2011 Gephi Consortium)
* [sigma.js](https://github.com/jacomyal/sigma.js/tree/master/plugins/sigma.layout.forceAtlas2), [Graphology](https://github.com/graphology/graphology-layout-forceatlas2/blob/master/iterate.js) (MIT, Guillaume Plique)
* [Anim-Wotmap](https://git.42l.fr/HugoTrentesaux/animwotmap) (Hugo Trentesaux)

The ForceAtlas2 paper was released under CC BY, Copyright 2014 Jacomy et al.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, version 3 of the License.  
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.  
You should have received a copy of the GNU Affero General Public License along with this program. If not, see https://www.gnu.org/licenses/.
