# ForceAtlas2 Rust Minimal

A fork of [ForceAtlas2 Rust](https://framagit.org/ZettaScript/forceatlas2-rs)
that compiles on stable Rust.  Supports *only* `f32` type coordinates, and has
dropped all optimizations dependent on nightly Rust.

## Known Issues

These issues were discovered at time of fork, and have yet to be remedied.

* The unit tests do not all pass.
* Not all combinations of features compile cleanly.
* Examples depend on system packages.

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
