//import Complex from './complex';
//import Filter from './filter';
//
//import {
//  runMultiFilter,
//  runMultiFilterReverse,
//  complex,
//  evaluatePhase,
//
//} from './utils';

import 'dart:math';

import 'complex.dart';
import 'filter.dart';
import 'utils.dart';

/**
 * Fir filter
 */
class FirFilter implements Filter {
  List<double> f;
  List<Complex> b;
  BufPoint z;
  FirFilter(List<double> filter)
      : f = filter,
        z = FirFilter.initZero(filter.length - 1),
        // note: coefficients are equal to input response
        b = List.generate(
            filter.length, (index) => Complex(filter[index], 0.0)) {}

  static BufPoint initZero(double cnt) {
    List<double> r = [];
    var i;
    for (i = 0; i < cnt; i++) {
      r.add(0);
    }
    return BufPoint(
      buf: r,
      pointer: 0,
    );
  }

  double doStep(double input, BufPoint d) {
    d.buf[d.pointer] = input;
    double out = 0;
    for (var cnt = 0; cnt < d.buf.length; cnt++) {
      out += (this.f[cnt] * d.buf[(d.pointer + cnt) % d.buf.length]);
    }
    d.pointer = (d.pointer + 1) % (d.buf.length);
    return out;
  }

  calcInputResponse(List<double> input) {
    var tempF = FirFilter.initZero(this.f.length - 1);
    return runMultiFilter(input, tempF,
        (double input, dynamic coeffs) => this.doStep(input, coeffs));
  }

  calcResponse({required double Fs, required double Fr}) {
    // z = exp(j*omega*pi) = cos(omega*pi) + j*sin(omega*pi)
    // z^-1 = exp(-j*omega*pi)
    // omega is between 0 and 1. 1 is the Nyquist frequency.
    var theta = -pi * (Fr / Fs) * 2;
    var h = Complex(
      0,
      0,
    );
    for (var i = 0; i < this.f.length - 1; i++) {
      h = h.add(this.b[i].mul(Complex(
            cos(theta * i),
            sin(theta * i),
          )));
    }
    var m = h.magnitude();
    var res = FilterResponse(
      magnitude: m,
      phase: h.phase(),
      dBmagnitude: 20 * log(m) * log10e,
    );
    return res;
  }

  responsePoint({required double Fs, required double Fr}) {
    return this.calcResponse(Fs: Fs, Fr: Fr);
  }

  List<FilterResponse> response(int? resolution) {
    resolution = resolution ?? 100;
    var r = resolution * 2;
    List<FilterResponse> res = List.generate(
        resolution,
        (index) => this.calcResponse(
              Fs: r.toDouble(), // TODO: Correct or should it be an int?
              Fr: index.toDouble(), // TODO: Correct or should it be an int?
            ));
    evaluatePhase(res);
    return res;
  }

  List<double> simulate(List<double> simInput) {
    return this.calcInputResponse(simInput);
  }

  double singleStep(double input) {
    return this.doStep(input, this.z);
  }

  List<double> multiStep(List<double> simInput, {bool overwrite = false}) {
    return runMultiFilter<BufPoint>(
        simInput, this.z, (input, coeffs) => this.doStep(input, coeffs),
        overwrite: overwrite);
  }

  filtfilt(dynamic input, {bool overwrite = false}) {
    return runMultiFilterReverse<BufPoint>(
        runMultiFilter<BufPoint>(
            input, this.z, (input, coeffs) => this.doStep(input, coeffs),
            overwrite: overwrite),
        this.z,
        (input, coeffs) => this.doStep(input, coeffs),
        overwrite: true);
  }

  reinit() {
    this.z = FirFilter.initZero(this.f.length - 1);
  }
}

class FilterResponse {
  double? magnitude;
  double? dBmagnitude;
  double? phase;
  double? unwrappedPhase;
  double? phaseDelay;
  double? groupDelay;
  FilterResponse({
    this.magnitude,
    this.dBmagnitude,
    this.phase,
    this.unwrappedPhase,
    this.phaseDelay,
    this.groupDelay,
  });
}

class BufPoint {
  List<double> buf;
  int pointer;
  BufPoint({required this.buf, required this.pointer});
}
