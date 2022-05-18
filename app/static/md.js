var chemobjs=[], nCanvas=3;

function resize_canvas() {
    var elem = $("#canvas-container-0"),
        width = elem.width(),
        height = Math.min(600, width/1.33333);
    
        chemobjs.forEach(function(chemobj) {
        chemobj.resize(width, height);
    });

}

function init_canvas() {
    for(var i=0; i<nCanvas; i++) {
        var chemobj = new ChemDoodle.TransformCanvas3D('canvas-viewer-'+i);
        chemobj.specs.set3DRepresentation('Stick');
        chemobj.specs.bonds_cylinderDiameter_3D = 3;
        chemobj.specs.atoms_sphereDiameter_3D = 3;
        chemobj.specs.backgroundColor = '#ffffff';
        chemobjs.push(chemobj);
    }
    resize_canvas();

    $.get(
        "static/examples.json",
        function(data) {
            load_file(data);
        }
    );
}

function load_file( data ) {
    for(var i=0; i<nCanvas; i++) {
        var molecule = ChemDoodle.readMOL(data[i]['mol'][0]);
        chemobjs[i].loadMolecule(molecule);
        $("#caption-"+(i+1)).text(data[i]['name'])
    }
}


$( document ).ready(function() {

    //-----------
    $(window).resize(function() {
        resize_canvas();
    });

    $("#copy-btn").click(function() {
        $.get(
            "static/examples.json",
            function(data) {
                current_mol = data[$('div.active').index()]
                $("#ror").val(current_mol.ror)
                $("#coo").val(current_mol.coo)
                $("#co").val(current_mol.co)
                $("#ch").val(current_mol.ch)
                $("#x").val(current_mol.x)
                $("#y").val(current_mol.y)
                $("#z").val(current_mol.z)
                $("#pf").val(current_mol.pf)
            }
        );
    })

    init_canvas();


    
});


