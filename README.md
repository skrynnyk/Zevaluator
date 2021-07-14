# Z-valuator

`Z-val` is a single-header C++11 library for efficient transformation and real-time evaluation of Discrete Transfer Functions (DTFs) in  the difference-equation form -- you just have to provide the numerator and denominator coefficients!

# Table of Contents

- Overview
	- Making a Discrete TF
	- Using a Discrete TF
	- Algebraic Operator/Utility DTFs
	- Algebraic Operator/Utility DTFs -- Integral / Derivative / Delay
	- Algebraic Operator/Utility DTFs -- Low-Pass Filters
		- Moving Average (FIR)
		- Exponential Smoothing (IIR)
	- Points of Configuration
- API Reference
- Benchmarks
- License


# Overview

Transfer functions are an incredibly flexible tool for describing LTI systems that are at the core of DSP and Control Theory domains.

At the heart, `Zval` is a small, efficient utility library that serves a very narrow purpose -- handling the boilerplate of implementing and evaluating a Discrete Transfer Function -- that specifically has resource-limited embedded targets in mind. 

For those just starting out in DSP or Control Theory, `Zval` also provides a number of helper utilities for creating moving average FIR and IIR low pass filters, derivative and integral operators, all requiring no deep mathematical knowledge.

### Making a Discrete TF
This is the general form you follow to create any DTF that you desire:
```c++
#include <zval.hpp>

auto Gz = zval::makeDTF(
    zval::makeCs(1.0f),         // Numerator coefs
    zval::makeCs(1.0f, -0.5f),  // Denominator coefs
    TickPolicy{},               // System tick getter
    100                         // Sampling time (in ticks)
);
```
You simply specify coefficients just as they appear in the discretized TF, along with sampling information. The only point worth noting is the `TickPolicy` class that must be defined so that the created `zval::DTF` object knows when to compute the next output value in the interval defined by the sampling time argument. For illustration, implementing a `TickPolicy` class for an Arduino project that wishes to use a time base in milliseconds, would look like this:
```c++
struct TickPolicy {
    static uint32_t getTick() { return millis(); }
};
```

### Using a Discrete TF
Using the resulting `zval::DTF` object `Gz` is very straightforward -- it's a function object, so you can just call it as if it were a function:
```c++
float input = 42;
float output = Gz(input);
```
> **Detail Warning:** the default interpolation behavior between sampling periods is the Zero-Order Hold (ZOH), meaning `Gz` output will remain at its previous value until the subsequent sampling period has elapsed.

### Algebraic Operator/Utility DTFs
There are a few special numerator and denominator coefficient sets for discrete transfer functions that create some general mathematical operators, such as: `derivative`, `integral`, `delay`, `moving average`, etc. This library provides a convenient way to create these operators and two very basic low pass filters through a set of utility functions:

- `dtfe::util::makeDerivative()`
- `dtfe::util::makeIntegral()`
-  `dtfe::util::makeDelay()`
-  `dtfe::util::makeMovingAvg()`
-  `dtfe::util::makeSinglePoleIIR()`

### Algebraic Operator/Utility DTFs -- Integral / Derivative / Delay
For example, if you need to integrate a signal every 10 ms using the Trapezoidal Rule, then you can simply create a trapezoidal integrator DTF and pass it your signal, like so:
```c++
// Size of `Gz_integrator` is only 32 bytes (using `float` for calculations)
auto Gz_integrator = zval::util::makeTrapezoidalIntegrator(TickPolicy{}, 10);

// Just emulating some data acquisition loop
while (is_reading) {
	Gz_integrator(readSignal());
}

// Parameter-less call returns the last computed DTF value
float signal_sum = Gz_integrator(); 
```
There are three possible methods of integration that you can represent with a DTF, and for which `Zval` provides helper constructor functions for:

 - `dtfe::util::makeForwardIntegrator()`
	 - "Forward Euler" or "Forward Rectangular" or "Left-hand Rule"
 - `dtfe::util::makeBackwardIntegrator()`
	 - "Backward Euler" or "Backward Rectangular" or "Right-hand Rule"
 - `dtfe::util::makeTrapezoidalIntegrator()`
	 - "Trapezoidal Rule"
	
> **Detail Warning:** integration of most signals tends to be unstable (i.e. grows unboundedly), so it's usually a good idea to clamp/saturate the integrator output at some applicable min/max value.

### Algebraic Operator/Utility DTFs -- Low-Pass Filters
If you are interested in simply smoothing out some of the high-frequency noise out of your sensor reading, and are not designing a bespoke filter, then you're likely looking for these two simple Low-Pass Filters (LPFs): 

 - Moving Average (FIR)
 - Exponential Smoothing (IIR)

These two filters are both very effective for simple "signal smoothing" use cases, but they both have relative Pros and Cons, briefly enumerated below:

| Filter | Stability | Phase Distortion | Roundoff Error |    Delay   | Memory Use | Compute Time |
|:------:|:---------:|:----------------:|:--------------:|:----------:|:----------:|:------------:|
|   FIR  |   Always  |      Uniform     |    Constant    | Med - High | Med - High |  Med - High  |
|   IIR  | Sometimes |      Varying     |   Accumulates  |     Low    |     Low    |      Low     |


#### | Moving Average (FIR)
This filter is...

#### | Exponential Smoothing (IIR)
This filter is...

### Points of Configuration
`Zval` is implemented with configurability in mind by using "policy" classes passed as template parameters. To go in hand with this, `Zval` provides a number of various policies for coefficient sets, platform-specific system ticks, and for DTF interpolation methods.

 

### Advanced Features
Should add:
 - Compile-time coefficient expansion for `N`-order delays
 - DTF wrappers, e.g. for implementation of true derivative

# Benchmarking

Some space+time data 

# API Reference

## `zval::makeCs()` 

```c++
template<class... TCs>
constexpr auto makeCs(const TCs&... cs) -> DistinctCSet<...>
```
This function does this and that

## `zval::makeUniformCs()` 

```c++
template<class... TCs>
constexpr auto makeUniformCs(const TCs&... cs) -> UniformCSet<...>
```
This function does this and that

## `zval::makeDTF()` 

```c++
template<class TAs,  class TBs,  class TickPolicy>
constexpr auto makeDTF(const TAs& As,  const TBs& Bs, TickPolicy, uint16_t Ts)
    -> DTF<...>
```
This function does this and that

## `zval::util::makeMovingAvg()` 

```c++
template<size_t TNum,  class TickPolicy>
constexpr auto makeMovingAvg(TickPolicy p, uint16_t Ts) -> DTF<...>
```
This function does this and that

## `zval::util::makeSinglePoleIIR()` 

```c++
template<class TickPolicy>
constexpr auto makeSinglePoleIIR(float alpha, TickPolicy p, uint16_t Ts) -> DTF<...>
```
This function does this and that

## `zval::util::makeDerivative()` 

```c++
template<class TickPolicy>
constexpr auto makeDerivative(TTickPolicy p, uint16_t Ts) -> DTF<...>
```
This function does this and that

## `zval::util::makeIntegral()` 

```c++
template<class TickPolicy>
constexpr auto makeIntegral(TTickPolicy p, uint16_t Ts) -> DTF<...>
```
This function does this and that

## `zval::util::makeDelay()` 

```c++
template<size_t TNum, class TickPolicy>
constexpr auto makeDelay(TTickPolicy p, uint16_t Ts) -> DTF<...>
```
This function does this and that

## License

Please see the `LICENSE` file for details.