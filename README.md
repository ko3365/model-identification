# Model (System) Identification
Consider a MIMO system - _m_ outputs, _q_ inputs, _ns_ states:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;x_{k&plus;1}&space;=&space;Ax_k&plus;Bu_k,&space;\quad&space;y_k&space;=&space;Cx_k}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;A\in&space;\mathbb{R}^{n_s\times&space;n_s}\quad&space;B\in&space;\mathbb{R}^{n_s\times&space;q}\quad&space;C\in&space;\mathbb{R}^{m\times&space;n_s}}">

### Impulse Response Data
Data provided is an input-output data provided by two-input/two-output system. u1_impulse.mat and u2_impulse.mat data with only one input channel is "on" in each data set. The state dimension is not known,
with the objective being the determination of a state-space discrete time system that can replicate the observed data.

```Matlab
load data/u1_impulse.mat
load data/u2_impulse.mat
% input channel 1 'on'
y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data;
% input channel 2 'on'
y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;
```
