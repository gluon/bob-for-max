# bob~ (Max external)

Moog ladder resonant filter simulation for Max/MSP.

Port of the **bob~** object from **Pure Data**, originally written by  **Miller Puckette**.
Platforms : macOS (intel, silicon) + win

The filter is implemented using **Runge-Kutta numerical integration** of a nonlinear differential equation modelling the Moog transistor ladder.

(thanks to Mathieu Staub for pointing me to bob~ during this class discussion)

References:

- Tim Stilson — *Analyzing the Moog VCF with Considerations for Digital Implementation*
- Timothy E. Stinchcombe — *Analysis of the Moog Transistor Ladder*
- Antti Huovilainen — *Non-linear digital implementation of the Moog ladder filter*

---

## Installation

1. Download the latest release: [bob-for-max-v1.0.0.zip](https://github.com/gluon/bob-for-max/releases)
2. Unzip the archive.
3. Place the `bob~` folder into your Max Packages directory:
   - **macOS**: `~/Documents/Max 9/Packages/`
   - **Windows**: `Documents\Max 9\Packages\`
4. Restart Max.

---

## Inlets

1. **signal** — audio input
2. **signal / float** — cutoff frequency (Hz)
3. **signal / float** — resonance

## Outlet

1. **signal** — filtered output

---

## Attributes

`@oversample`  
Oversampling factor for the Runge-Kutta solver.

`@saturation`  
Controls the transistor saturation function.

---

## Messages

`clear`  
Reset filter state.

---

## License

This software is based on the **bob~ object from Pure Data** by **Miller Puckette** and others.

Released under the **Standard Improved BSD License**.  
See `LICENSE.txt`.

---

## Author

Max port by **Julien Bayle** / Structure Void

https://structure-void.com
