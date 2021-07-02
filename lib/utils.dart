import 'dart:math';

import 'complex.dart';

/**
 * Evaluate phase
 */

evaluatePhase(List<dynamic> res) {
  var tpi = 2 * pi;
  List<double> phase = [];
  for (var cnt = 0; cnt < res.length; cnt++) {
    phase.add(res[cnt].phase);
  }
  res[0].unwrappedPhase = res[0].phase;
  res[0].groupDelay = 0.0; // TODO: was 0
  // TODO: more sophisticated phase unwrapping needed
  for (var cnt = 1; cnt < phase.length; cnt++) {
    var diff = phase[cnt] - phase[cnt - 1];
    if (diff > pi) {
      for (var xcnt = cnt; xcnt < phase.length; xcnt++) {
        phase[xcnt] -= tpi;
      }
    } else if (diff < -pi) {
      for (var xcnt = cnt; xcnt < phase.length; xcnt++) {
        phase[xcnt] += tpi;
      }
    }
    if (phase[cnt] < 0) {
      res[cnt].unwrappedPhase = -phase[cnt];
    } else {
      res[cnt].unwrappedPhase = phase[cnt];
    }

    res[cnt].phaseDelay = res[cnt].unwrappedPhase / (cnt / res.length);
    res[cnt].groupDelay =
        (res[cnt].unwrappedPhase - res[cnt - 1].unwrappedPhase) /
            (pi / res.length);
    if (res[cnt].groupDelay < 0) {
      res[cnt].groupDelay = -res[cnt].groupDelay;
    }
  }
  if (res[0].magnitude != 0) {
    res[0].phaseDelay = res[1].phaseDelay;
    res[0].groupDelay = res[1].groupDelay;
  } else {
    res[0].phaseDelay = res[2].phaseDelay;
    res[0].groupDelay = res[2].groupDelay;
    res[1].phaseDelay = res[2].phaseDelay;
    res[1].groupDelay = res[2].groupDelay;
  }
}

/**
 * Run multi filter
 */

List<double> runMultiFilter<T>(
    List<double> input, T d, double Function(double input, T d) doStep,
    {bool overwrite = false}) {
  List<double> out;
  if (overwrite) {
    out = input;
  } else {
    out = List.filled(input.length, 0.0);
  }
  for (var i = 0; i < input.length; i++) {
    out[i] = doStep(input[i], d);
  }
  return out;
}

List<double> runMultiFilterReverse<T>(
    List<double> input, T d, double Function(double input, T d) doStep,
    {bool overwrite = false}) {
  List<double> out;
  if (overwrite) {
    out = input;
  } else {
    out = List.filled(input.length, 0.0);
  }
  for (var i = input.length - 1; i >= 0; i--) {
    out[i] = doStep(input[i], d);
  }
  return out;
}

int factorial(int n, {int? a}) {
  a = a ?? 1;
  if (n != n.floor() || a != a.floor()) {
    return 1;
  }
  if (n == 0 || n == 1) {
    return a;
  } else {
    return factorial(n - 1, a: a * n);
  }
}

/**
 * Bessel factors
 */

List<int> besselFactors(int n) {
  List<int> res = [];
  for (var k = 0; k < n + 1; k++) {
    var p = factorial(2 * n - k);

    var q = pow(2, n - k) * factorial(k) * factorial(n - k);
    res.insert(0, (p / q).floor());
  }
  return res;
}
