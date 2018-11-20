"use strict";

import 'select2/dist/css/select2.css';
import 'select2-theme-bootstrap4/dist/select2-bootstrap.css';
import 'nouislider/distribute/nouislider.min.css';
import './css/style.scss';
import './css/nv.d3.css';

import 'bootstrap';
import 'select2/dist/js/select2.full.js';
import './js/my-fontawesome.js';
import './js/nv.d3.js';

import './js/project.js';
import './js/protein_page.js';
import './js/datatables.js';
import './js/favit.js';
import './js/jquery.formset.js';
import './js/onload.js';
import './js/microscope.js';
import LineageChart from './js/lineage.js'
window.LineageChart = LineageChart;

import initFRET from './js/fret.js';
window.initFRET = initFRET;

import initSpectra from './js/spectra/spectra.js';
window.initSpectra = initSpectra;

import FPPropChart from './js/ichart.js';
window.FPPropChart = FPPropChart;
