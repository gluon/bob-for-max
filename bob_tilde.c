#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
#include <math.h>

/*
bob~ — Moog analog lowpass resonant filter simulation

ported to Max by Julien Bayle from the Pure Data object bob~
original design and implementation: Miller Puckette

Imitates a Moog resonant filter by Runge-Kutta numerical integration
of a differential equation describing the transistor ladder dynamics.

References:

Tim Stilson — Analyzing the Moog VCF with Considerations for Digital Implementation
https://ccrma.stanford.edu/~stilti/papers/moogvcf.ps.gz

Timothy E. Stinchcombe — Analysis of the Moog Transistor Ladder and Derivative Filters

Antti Huovilainen — Non-linear digital implementation of the moog ladder filter

Differential equations:

y1' = k * (S(x - r * y4) - S(y1))
y2' = k * (S(y1) - S(y2))
y3' = k * (S(y2) - S(y3))
y4' = k * (S(y3) - S(y4))

k = cutoff frequency coefficient
r = resonance feedback (<= 4 stable)
S(x) = saturation function

Inlets:
1. signal input
2. cutoff frequency (Hz) — signal or float
3. resonance — signal or float

Outlet:
1. filtered signal

License:
Standard Improved BSD License (see LICENSE.txt)
*/

#define DIM 4
#define FLOAT double

typedef struct _params {
    FLOAT p_input;
    FLOAT p_cutoff;
    FLOAT p_resonance;
    FLOAT p_saturation;
} t_params;

typedef struct _bob {
    t_pxobject x_obj;
    t_params   x_params;
    FLOAT      x_state[DIM];
    double     x_sr;
    int        x_oversample;
    double     default_cutoff;
    double     default_resonance;
    short      x_cutoff_connected;
    short      x_res_connected;
} t_bob;

static t_class *bob_class;

void bob_assist(t_bob *x, void *b, long m, long a, char *s);

static inline FLOAT clip(FLOAT value, FLOAT saturation, FLOAT satinv) {
    FLOAT v2 = (value * satinv > 1.0 ? 1.0 : (value * satinv < -1.0 ? -1.0 : value * satinv));
    return (saturation * (v2 - (1.0/3.0) * v2 * v2 * v2));
}

static void calc_derivatives(FLOAT *dstate, FLOAT *state, t_params *params) {

    FLOAT k = (2.0 * M_PI) * params->p_cutoff;

    FLOAT sat = params->p_saturation;
    FLOAT satinv = 1.0 / sat;

    FLOAT satstate0 = clip(state[0], sat, satinv);
    FLOAT satstate1 = clip(state[1], sat, satinv);
    FLOAT satstate2 = clip(state[2], sat, satinv);

    dstate[0] = k * (clip(params->p_input - params->p_resonance * state[3], sat, satinv) - satstate0);
    dstate[1] = k * (satstate0 - satstate1);
    dstate[2] = k * (satstate1 - satstate2);
    dstate[3] = k * (satstate2 - clip(state[3], sat, satinv));
}

static void solver_rungekutte(FLOAT *state, FLOAT stepsize, t_params *params) {

    int i;

    FLOAT deriv1[DIM];
    FLOAT deriv2[DIM];
    FLOAT deriv3[DIM];
    FLOAT deriv4[DIM];
    FLOAT tempstate[DIM];

    calc_derivatives(deriv1, state, params);

    for (i = 0; i < DIM; i++)
        tempstate[i] = state[i] + 0.5 * stepsize * deriv1[i];

    calc_derivatives(deriv2, tempstate, params);

    for (i = 0; i < DIM; i++)
        tempstate[i] = state[i] + 0.5 * stepsize * deriv2[i];

    calc_derivatives(deriv3, tempstate, params);

    for (i = 0; i < DIM; i++)
        tempstate[i] = state[i] + stepsize * deriv3[i];

    calc_derivatives(deriv4, tempstate, params);

    for (i = 0; i < DIM; i++)
        state[i] += (1.0/6.0) * stepsize * (deriv1[i] + 2.0*deriv2[i] + 2.0*deriv3[i] + deriv4[i]);
}

void bob_perform64(
    t_bob *x,
    t_object *dsp64,
    double **ins,
    long numins,
    double **outs,
    long numouts,
    long sampleframes,
    long flags,
    void *userparam)
{

    double *in         = ins[0];
    double *cutoff_sig = ins[1];
    double *res_sig    = ins[2];

    double *out = outs[0];

    FLOAT s0 = x->x_state[0];
    FLOAT s1 = x->x_state[1];
    FLOAT s2 = x->x_state[2];
    FLOAT s3 = x->x_state[3];

    FLOAT stepsize = 1.0 / (x->x_oversample * x->x_sr);

    t_params p = x->x_params;

    for (int i = 0; i < sampleframes; i++) {

        p.p_input = in[i];

        if (x->x_cutoff_connected)
            p.p_cutoff = (cutoff_sig[i] < 20.0) ? 20.0 : cutoff_sig[i];
        else
            p.p_cutoff = x->default_cutoff;

        if (x->x_res_connected)
            p.p_resonance = (res_sig[i] < 0.0) ? 0.0 : res_sig[i];
        else
            p.p_resonance = x->default_resonance;

        FLOAT current_state[4] = {s0, s1, s2, s3};

        for (int j = 0; j < x->x_oversample; j++) {
            solver_rungekutte(current_state, stepsize, &p);
        }

        s0 = current_state[0];
        s1 = current_state[1];
        s2 = current_state[2];
        s3 = current_state[3];

        if (isnan(s3) || isinf(s3)) {
            s0 = 0.0;
            s1 = 0.0;
            s2 = 0.0;
            s3 = 0.0;
        }

        out[i] = s3;
    }

    x->x_state[0] = s0;
    x->x_state[1] = s1;
    x->x_state[2] = s2;
    x->x_state[3] = s3;
}

void bob_dsp64(
    t_bob *x,
    t_object *dsp64,
    short *count,
    double samplerate,
    long maxvectorsize,
    long flags)
{

    x->x_sr = samplerate;

    x->x_cutoff_connected = count[1];
    x->x_res_connected    = count[2];

    dsp_add64(
        dsp64,
        (t_object *)x,
        (t_perfroutine64)bob_perform64,
        0,
        NULL
    );
}

void bob_float(t_bob *x, double f) {

    int inlet = proxy_getinlet((t_object *)x);

    if (inlet == 1)
        x->default_cutoff = f;

    else if (inlet == 2)
        x->default_resonance = f;
}

void bob_clear(t_bob *x) {

    for (int i = 0; i < DIM; i++)
        x->x_state[i] = 0.0;
}

void *bob_new(t_symbol *s, long argc, t_atom *argv) {

    t_bob *x = (t_bob *)object_alloc(bob_class);

    if (x) {

        dsp_setup((t_pxobject *)x, 3);

        outlet_new((t_object *)x, "signal");

        x->x_sr = sys_getsr();

        if (x->x_sr <= 0)
            x->x_sr = 44100.0;

        x->x_oversample = 2;

        x->x_params.p_saturation = 3.0;

        x->default_cutoff =
            (argc > 0) ? atom_getfloat(&argv[0]) : 500.0;

        x->default_resonance =
            (argc > 1) ? atom_getfloat(&argv[1]) : 0.5;

        x->x_cutoff_connected = 0;
        x->x_res_connected = 0;

        bob_clear(x);
    }

    return x;
}

void bob_free(t_bob *x) {

    dsp_free((t_pxobject *)x);
}

void bob_assist(t_bob *x, void *b, long m, long a, char *s)
{
    if (m == ASSIST_INLET) {
        switch (a) {
            case 0: snprintf(s, 256, "(signal) input"); break;
            case 1: snprintf(s, 256, "(signal/float) cutoff frequency (Hz)"); break;
            case 2: snprintf(s, 256, "(signal/float) resonance"); break;
        }
    } else {
        snprintf(s, 256, "(signal) filtered output");
    }
}

C74_EXPORT void ext_main(void *r) {

    t_class *c =
        class_new(
            "bob~",
            (method)bob_new,
            (method)bob_free,
            (long)sizeof(t_bob),
            0L,
            A_GIMME,
            0
        );

    class_addmethod(c,
        (method)bob_dsp64,
        "dsp64",
        A_CANT,
        0);

    class_addmethod(c,
        (method)bob_float,
        "float",
        A_FLOAT,
        0);

    class_addmethod(c,
        (method)bob_clear,
        "clear",
        0);

    class_addmethod(c,
        (method)bob_assist,
        "assist",
        A_CANT,
        0);

    CLASS_ATTR_LONG(c,
        "oversample",
        0,
        t_bob,
        x_oversample);

    CLASS_ATTR_FILTER_MIN(c,
        "oversample",
        1);

    CLASS_ATTR_DOUBLE(c,
        "saturation",
        0,
        t_bob,
        x_params.p_saturation);

    CLASS_ATTR_FILTER_MIN(c,
        "saturation",
        0.001);

    class_dspinit(c);

    class_register(CLASS_BOX, c);

    bob_class = c;
}
