import 'smiles-drawer/dist/smiles-drawer.min.js';
import LiteMol from './js/LiteMol-Plugin.js';
import './css/litemol/LiteMol-plugin-blue.css'


var smilesDrawer = new SmilesDrawer.Drawer({
    compactDrawing: false,
    height: 300,
});

const PluginSpec = LiteMol.Plugin.getDefaultSpecification();
const LayoutRegion = LiteMol.Bootstrap.Components.LayoutRegion;
const Components = LiteMol.Plugin.Components;
PluginSpec.components = [
  Components.Visualization.HighlightInfo(LayoutRegion.Main, true),
  Components.Entity.Current('LiteMol', LiteMol.Plugin.VERSION.number)(LayoutRegion.Right, true),
  Components.Transform.View(LayoutRegion.Right),
  //Components.Context.Log(LayoutRegion.Bottom, true),
  Components.Context.Overlay(LayoutRegion.Root),
  //Components.Context.Toast(LayoutRegion.Main, true),
  Components.Context.BackgroundTasks(LayoutRegion.Main, true)
];

function initLiteMol(selection){

    var plugin = LiteMol.Plugin.create({
        customSpecification: PluginSpec,
        target: selection,
        viewportBackground: '#fff',
        layoutState: {
            hideControls: true,
            isExpanded: false
        },
        allowAnalytics: true
    });
    // To see all the available methods on the SimpleController,
    // please check src/Plugin/Plugin/SimpleController.ts
    //////////////////////////////////////////////////////////////
    //
    // The underlaying instance of the plugin can be accessed by
    //
    //   plugin.instance
    //////////////////////////////////////////////////////////////
    //
    // To create and apply transforms, use
    //
    //   let t = plugin.createTransform();
    //   t.add(...).then(...);
    //   plugin.applyTransform(t);
    //
    // Creation of transforms is illusted in other examples.
    //////////////////////////////////////////////////////////////
    //
    // To execute commands, the SimpleController provides the method command.
    //
    //   plugin.command(command, params);
    //
    // To find examples of commands, please see the Commands example.
    //////////////////////////////////////////////////////////////
    //
    // To subscribe for events, the SimpleController provides the method subscribe.
    //
    //   plugin.subscribe(event, callback);
    //
    // To find examples of events, please see the Commands example as well.
    // It shows how to subscribe interaction events, where available events are located, etc.

    return plugin
}


function getPDBbinary (id) {
  return new Promise(function (resolve, reject) {
    var xhr = new XMLHttpRequest();
    xhr.open("GET", "https://www.ebi.ac.uk/pdbe/static/entry/" + id.toLowerCase() + "_updated.cif");
    xhr.onload = function () {
      if (this.status >= 200 && this.status < 300) {
        resolve(xhr.response);
      } else {
        reject({
          status: this.status,
          statusText: xhr.statusText
        });
      }
    };
    xhr.onerror = function () {
      reject({
        status: this.status,
        statusText: xhr.statusText
      });
    };
    xhr.send();
  });
}


function downloadPDBMeta(pdbidstring){
  return $.when(
    $.get({
        url: "https://www.rcsb.org/pdb/rest/describeMol?structureId=" + pdbidstring,
        dataType: "xml",
        success: function(xml){
          $(xml).find('structureId').each(function(i, p){
            var el = $(p);
            var pdbid = el.attr('id');
            //pdb_info[pdbid].macroMolecule = el.find('macroMolecule').attr('name');
            pdb_info[pdbid].chains = el.find('chain').map(function(i, d){return d.id}).toArray();
          })
        }
    }),
    $.get({
        url: "https://www.rcsb.org/pdb/rest/describePDB?structureId=" + pdbidstring,
        dataType: "xml",
        success: function(xml){
          $(xml).find('PDB').each(function(i, p){
            var pdb = $(p);
            var pdbid = pdb.attr('structureId');
            pdb_info[pdbid].title = pdb.attr('title');
            pdb_info[pdbid].expMethod = pdb.attr('expMethod');
            pdb_info[pdbid].resolution = parseFloat(pdb.attr('resolution'));
            pdb_info[pdbid].deposition_date = pdb.attr('deposition_date');
            pdb_info[pdbid].last_modification_date = pdb.attr('last_modification_date');
            pdb_info[pdbid].structure_authors = pdb.attr('structure_authors');
            pdb_info[pdbid].pubmedId = pdb.attr('pubmedId');
          })
        }
    }),
    $.get({
        url: "https://www.rcsb.org/pdb/rest/ligandInfo?structureId=" + pdbidstring,
        dataType: "xml",
        success: function(xml){
          $(xml).find('structureId').each(function(i, p){
            var el = $(p);
            var pdbid = el.attr('id');
            var heavy = el.find('ligand').toArray().reduce(function (a, b) {
              return +$(a).attr('molecularWeight') > +$(b).attr('molecularWeight') ? a : b;
            });
            pdb_info[pdbid].smiles = $(heavy).find('smiles').text().replace(/\\\\/g, '\\').replace('CC(=O)O', '*');
            pdb_info[pdbid].chemId = $(heavy).attr('chemicalID');
            pdb_info[pdbid].chemicalName = $(heavy).find('chemicalName').text();
            pdb_info[pdbid].formula = $(heavy).find('formula').text();
            pdb_info[pdbid].molecularWeight = $(heavy).attr('molecularWeight');

          })
        }
    })
    )
}

function loadSmiles(pdbid){
  var canv = $('#smilesCanvas');
  canv[0].getContext('2d').clearRect(0, 0, canv.width, canv.height);
  if (canv.data('chemId') != pdb_info[pdbid].chemId){
    SmilesDrawer.parse(pdb_info[pdbid].smiles, function(tree) {
      canv.remove();
      canv = $('<canvas>', {id: 'smilesCanvas'}).appendTo($('#smilesDrawing .card-img-top'));
      smilesDrawer.draw(tree, 'smilesCanvas', 'light', false);
      canv.data('chemId', pdb_info[pdbid].chemId)
      canv.height('auto');
      canv.width('100%');
    });
  }
}

function loadChemInfo(pdbid){
  var e = pdb_info[pdbid];
  $('#chem-title').html(e.title);
  $('#chem-date').html(e.deposition_date);
  $('#chem-authors').html(e.structure_authors.split(',')[0] + ' et al. ');
  $('#chem-pubmed').attr('href', 'https://www.ncbi.nlm.nih.gov/pubmed/' + e.pubmedId);
  $('#chem-id').html('<a target="_blank" class="text-secondary" href="https://www.rcsb.org/ligand/' + e.chemId + '">' + e.chemId + '</a>')
  $('#chem-form').html(e.formula)
}

//variable to store LiteMol Component Scope which has all the methods
var liteMolScope;
var pdb_info = {}
var getPDBinfo = function(pdbids){
  if (Array.isArray(pdbids)) {
    for (var i = 0; i < pdbids.length; ++i) {
      pdbids[i] = pdbids[i].toUpperCase();
      pdb_info[pdbids[i]] = {}
    }
    var pdbidstring = pdbids.join(",");
  } else {
    pdbids = pdbids.toUpperCase();
    pdb_info[pdbids] = {};
    var pdbidstring = pdbids;
  }

  downloadPDBMeta(pdbidstring).done(function(){

      var select = $("#pdb_select");
      $(pdbids.sort(function(a, b){
        return pdb_info[a].resolution - pdb_info[b].resolution
      })).each(function(i, d){
        select.append($('<option>', {value: d}).html(
          `${d} (${pdb_info[d].resolution} Å, ${pdb_info[d].chains.length}-chain)`
        ))
      })

      var plugin = initLiteMol('#litemol-viewer')
      plugin.dataCache = {};

      select.change(function(){
        var id = this.value;
        if (!plugin.dataCache.hasOwnProperty(id)){ plugin.dataCache[id] = getPDBbinary(id); }
        plugin.clear();
        plugin.dataCache[id].then(data => plugin.loadMolecule({ data, id }));
        $("#pdb-out-link").attr('href', "https://www.rcsb.org/structure/" + this.value);
        if (pdb_info[this.value]){
          loadSmiles(this.value);
          loadChemInfo(this.value);
        }
      }).trigger('change');

      $('body').click(function(e) {
        if ($(".lm-layout-right").length){
          if ($(e.target).closest('#litemol-viewer').length === 0) {
              plugin.setLayoutState({hideControls: true});
          }
        }
      });
   })
   .fail(function(){
      $(pdbids).each(function(i, d){
        $("#pdb_select").append($('<option>', {value: d}).html(d))
      })
      $("#pdb-info").html('<small class="text-muted">failed to retrieve metadata from PDB</small>');
   });
}

var loadDensity = function() {
  liteMolScope.LiteMolComponent.loadDensity()
  $("#density_loader").text('toggle electron density')
  $("#density_loader").attr('onclick', 'ToggleDensity()')
}

var ToggleDensity = function() {
  liteMolScope.LiteMolComponent.toggleDensity()
}

export default function initPDB(pdbids) {
  getPDBinfo(pdbids);
}

window.initPDB = initPDB;
