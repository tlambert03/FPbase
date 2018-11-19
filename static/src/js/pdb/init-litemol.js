/*
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
var LiteMol;
(function (LiteMol) {
    var SimpleControllerExample;
    (function (SimpleControllerExample) {
        var plugin = LiteMol.Plugin.create({
            
            target: '#app',
            viewportBackground: '#fff',
            layoutState: {
                hideControls: true,
                isExpanded: false
            },
            allowAnalytics: true
        });
        var id = '1gfl';
        plugin.loadMolecule({
            id: id,
            format: 'cif',
            url: "https://www.ebi.ac.uk/pdbe/static/entry/" + id.toLowerCase() + "_updated.cif",
            // instead of url, it is possible to use
            // data: "string" or ArrayBuffer (for BinaryCIF)
            // loaded molecule and model can be accessed after load
            // using plugin.context.select(modelRef/moleculeRef)[0],
            // for example plugin.context.select('2yog-molecule')[0]
            moleculeRef: id + '-molecule',
            modelRef: id + '-model',
        }).then(function () {
            // Use this (or a modification of this) for custom visualization:
            // const style = LiteMol.Bootstrap.Visualization.Molecule.Default.ForType.get('BallsAndSticks');  
            // const t = plugin.createTransform();
            // t.add(id + '-model', LiteMol.Bootstrap.Entity.Transformer.Molecule.CreateVisual, { style: style })
            // plugin.applyTransform(t);
            console.log('Molecule loaded');
        }).catch(function (e) {
            console.error(e);
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
    })(SimpleControllerExample = LiteMol.SimpleControllerExample || (LiteMol.SimpleControllerExample = {}));
})(LiteMol || (LiteMol = {}));
