export const options = {
    showArea: true,
    minwave: 300,
    maxwave: 1000,
    startingBrush: [350, 800],
    autoscaleBrush: false,
    exNormWave: undefined,
    scale: 'linear',
    hide2p: true,
    scaleToEC: false,
    scaleToQY: false,
};

export const userOptions = {
    showArea: {type: 'checkbox', msg: 'Fill area under curve'},
    autoscaleBrush: {type: 'checkbox', msg: 'Auto-rescale X-axis (using zoom above auto-disables this)'},
    hide2p: {type: 'checkbox', msg: 'Hide 2-photon spectra by default'},
    scaleToEC: {type: 'checkbox', msg: 'Scale excitation spectra to extinction coefficient (% of highest fluor)'},
    scaleToQY: {type: 'checkbox', msg: 'Scale emission spectra to quantum yield'},
};

export const CONST = {
    category: {
        dye: 'd',
        protein: 'p',
        light: 'l',
        filter: 'f',
        camera: 'c',
    },
    stype: {
        ex: 'ex',
        abs: 'ab',
        em: 'em',
        twop: '2p',
        bp: 'bp',
        bpx: 'bx',
        bpm: 'bm',
        sp: 'sp',
        lp: 'lp',
        bs: 'bs',
        qe: 'qe',
        pd: 'pd',
    }
};
