

# Z-evaluator

`Z-eval` is a single-header C++11 library for efficient transformation and real-time evaluation of Discrete Transfer Functions (DTFs) in the difference-equation form -- you simply have to specify the numerator and denominator coefficients.

# Table of Contents

- Overview
    - Making a Discrete TF
    - Using a Discrete TF
    - Algebraic/Utility DTFs
      - Integral / Derivative / Delay
      - Low-Pass Filters
      - PID Controller
    - Points of Configuration
    - Advanced Features
- API Reference
- Benchmarks
- License


# Overview

Transfer functions are an incredibly helpful tool for describing Linear Time Invariant (LTI) systems that make up a key part of Control Theory and Digital Signals Processing.

At its core, `Z-eval` is a small and efficient utility serving a very narrow purpose -- abstracting away the boilerplate of implementing and evaluating difference equation forms of Discrete Transfer Functions. This utility specifically focuses on resource-limited embedded targets. 

For those just starting out in DSP or Control Theory, `Z-eval` also provides a few helper functions for creating Moving Average and IIR filters, derivative and integral approximations, and group delays without requiring knowledge of any of the underlying arithmetic. 

### Creating a Discrete TF
This is the general form you follow to create any DTF that you desire:
```c++
#include <zeval.hpp>

auto Gz = zeval::makeDTF(
    zeval::makeCs(1.0f),         // Numerator coefficients
    zeval::makeCs(1.0f, -0.5f),  // Denominator coefficients
    TickPolicy{},                // System tick policy object
    100                          // Sampling time (in TickPolicy tick units)
);
```
You simply specify coefficients as they appear in your discretized TF, along with the sampling time used for the discretization. The only point worth noting is the `TickPolicy` class that must be defined by the user such that the created `zeval::DTF` object knows when to compute the next output value in the interval defined by the sampling time argument.

As an example, implementing a `TickPolicy` class for an Arduino project that wishes to use a time base in milliseconds would simply look like this:
```c++
struct TickPolicy {
    static uint32_t getTick() const { return millis(); }
};
```

### Using a Discrete TF
Using the resulting `zeval::DTF` object `Gz` is very straightforward -- it's a function object, so you can simply call it as if it were a function:
```c++
float input = 42;
float output = Gz(input);
```
> **Detail:** the default interpolation behavior between sampling periods is the Zero-Order Hold (ZOH), meaning `Gz` output will remain at its previous value until the subsequent sampling period has elapsed. Specifying custom behavior is in scope for future development.

### Algebraic/Utility DTFs
There are a few special numerator and denominator coefficient sets for discrete transfer functions that create some general mathematical operators, such as: `derivative`, `integral`, `delay`, `moving average`, etc. This library provides a convenient set of utility functions to create these operators along with two very basic low pass filters:

- `dtfe::util::makeDerivative()`
- `dtfe::util::makeIntegral()`
- `dtfe::util::makeDelay()`
- `dtfe::util::makeMovingAvg()`
- `dtfe::util::makeSinglePoleIIR()`

#### Integral / Derivative / Delay
For example, if you need to integrate a signal every 10 ms using the Trapezoidal Rule, then you can simply create a trapezoidal integrator DTF and pass it your signal, like so:
```c++
// Size of `Gz_integrator` is only 32 bytes (using `float` for calculations)
auto Gz_integrator = zeval::util::makeTrapezoidalIntegrator(TickPolicy{}, 10);

// Emulating some data acquisition loop
while (is_reading) {
    Gz_integrator(readSignal());
}

// Tip: parameter-less call returns the last computed DTF value
float signal_sum = Gz_integrator(); 
```
There are several possible methods of integration that you can represent with a DTF, and for some of which `Z-eval` provides helper constructor functions for:

 - `zeval::util::makeForwardIntegrator()`
 - `zeval::util::makeBackwardIntegrator()`
 - `zeval::util::makeTrapezoidalIntegrator()`

#### Low-Pass Filters
If you want to smooth out your noisy sensor reading or some other signal, and do no need to or want to design a custom filter, then you're likely looking for one of these two very simple Low-Pass Filters (LPFs): 

 - Moving Average (FIR)
 - Exponential Smoothing (IIR)

These two LPF methods are effective for most simple "signal smoothing" use cases, but they both have their Pros and Cons, some of which are enumerated below:

| Filter |   Stability    | Phase Delay | Roundoff Error | Memory Use | Compute Time |
|:------:|:--------------:|:-----------:|:--------------:|:----------:|:------------:|
|   FIR  |   Guaranteed   |   Uniform   |    Constant    | Low - High |    Low - High    |
|   IIR  | Not guaranteed |   Varying   |   Accumulates  |     Lowest    |    Lowest     |

#### PID Controller
You can implement a PID controller in the form of a DTF, so `Z-eval` provides a helper function to do the necessary coefficient transformations for you so that you can focus on using it instead, which could look something like this:
```c++
// Transfer function of a tuned discrete-time PID controller
auto Gz_pid = zeval::util::makePID(
    1.9f,          // Proportional gain term
    0.5f,          // Integral gain term
    0.1f,          // Derivative gain term
    TickPolicy{},  // System tick policy object
    10             // Sampling time (in TickPolicy tick units)
);

float setpoint = 42.0f;

// Emulating the execution loop
while (is_running) {
    // Evaluate PID controller output given some sensor reading
    float error = setpoint - readSignal()
    float pid_output = Gz_pid(error);
    // Update the output of some actuator based on PID output
    setActuator(pid_output);
}
``` 
However, there are limitations for PID controllers implemented in a transfer function form. For example, there is no way to implement clamping integral anti-windup as it cannot be expressed as an LTI system -- clamping is a non-linearity.

### Points of Configuration
`Zeval` is implemented with configurability in mind by using "policy" classes passed as template parameters. `Zeval` provides a number of various policies for coefficient sets, platform-specific system ticks, and for DTF interpolation methods.

 

### Advanced Features
TODO:
 - Compile-time coefficient expansion for `N`-order delays
 - DTF wrappers, e.g. for implementation of true derivative

# Benchmarking

TODO: 
 - Space+time data 

# API Reference

## `zeval::makeCs()`

```c++
template<class... TCs>
constexpr auto makeCs(const TCs&... cs) -> DistinctCSet<...>
```
This function is used to produce a set of distinct coefficients, used for
specifying denominator and numerator TF coefficients.

## `zeval::makeUniformCs()`

```c++
template<class... TCs>
constexpr auto makeUniformCs(const TCs&... cs) -> UniformCSet<...>
```
This function is used to produce a set of coefficients of same value, used for
specifying denominator and numerator TF coefficients.

## `zeval::makeDTF()`

```c++
template<class TAs,  class TBs,  class TickPolicy>
constexpr auto makeDTF(const TAs& As,  const TBs& Bs, TickPolicy, uint16_t Ts)
    -> DTF<...>
```
This function creates a callable object, representing the DTF with provided 
coefficients, that is used to perform the difference equation evaluation

## `zeval::util::makeMovingAvg()`

```c++
template<size_t TNum,  class TickPolicy>
constexpr auto makeMovingAvg(TickPolicy p, uint16_t Ts) -> DTF<...>
```
This is a utility function for quickly generating a moving average (also 
sometimes called "boxcar" or "FIR") low pass filter DTF object.

## `zeval::util::makeSinglePoleIIR()`

```c++
template<class TickPolicy>
constexpr auto makeSinglePoleIIR(float alpha, TickPolicy p, uint16_t Ts)
    -> DTF<...>
```
This is a utility function for quickly generating a single-pole infinite 
impulse response (also sometimes called "recursive") low pass filter DTF object.

## `zeval::util::makeDerivative()`

```c++
template<class TickPolicy>
constexpr auto makeDerivative(TTickPolicy p, uint16_t Ts) -> DTF<...>
```
This is a utility function for quickly generating a derivative approximating 
DTF object.

## `zeval::util::makeIntegral()`

```c++
template<class TickPolicy>
constexpr auto makeIntegral(TTickPolicy p, uint16_t Ts) -> DTF<...>
```
This is a utility function for quickly generating an integral approximating 
DTF object.

## `zeval::util::makeDelay()`

```c++
template<size_t TNum, class TickPolicy>
constexpr auto makeDelay(TTickPolicy p, uint16_t Ts) -> DTF<...>
```
This is a utility function for quickly generating a signal-delaying DTF object.


## `zeval::util::makePID()`

```c++
template<class TTickPolicy>
constexpr auto makePID(float kp, float ki, float kd, TTickPolicy p, float Ts) 
    -> DTF<...>
```
This is a utility function for quickly generating a PID controller DTF object

# License

Please see the `LICENSE` file for details.