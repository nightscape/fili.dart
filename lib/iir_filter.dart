import 'dart:math';

import 'package:complex/complex.dart';
import 'package:fili.dart/complex_ext.dart';
import 'package:fili.dart/fir_filter.dart';

import 'filter.dart';
import 'iir_coeffs.dart';
import 'utils.dart';

// params: array of biquad coefficient objects and z registers
// stage structure e.g. {k:1, a:[1.1, -1.2], b:[0.3, -1.2, -0.4], z:[0, 0]}
class IirFilter implements Filter {
  List<TempF> cf;
  List<CC> cc;
  Complex cone = Complex(
    1,
    0,
  );
  late List<Coeffs> f;

  IirFilter(List<Coeffs> filter)
      : f = filter,
        cf = List.generate(
            filter.length,
            (cnt) => TempF(
                  b0: Complex(filter[cnt].b[0], 0),
                  b1: Complex(filter[cnt].b[1], 0),
                  b2: Complex(filter[cnt].b[2], 0),
                  a1: Complex(filter[cnt].a[0], 0),
                  a2: Complex(filter[cnt].a[1], 0),
                  k: Complex(filter[cnt].k, 0),
                  z: [0, 0],
                )),
        cc = List.generate(
            filter.length,
            (cnt) => CC(
                  b1: filter[cnt].b[1] / filter[cnt].b[0],
                  b2: filter[cnt].b[2] / filter[cnt].b[0],
                  a1: filter[cnt].a[0],
                  a2: filter[cnt].a[1],
                ));

  runStage(TempF s, double input) {
    var temp = input * s.k.real - s.a1.real * s.z[0] - s.a2.real * s.z[1];
    var out = s.b0.real * temp + s.b1.real * s.z[0] + s.b2.real * s.z[1];
    s.z[1] = s.z[0];
    s.z[0] = temp;
    return out;
  }

  double doStep(double input, List<TempF> coeffs) {
    var out = input;
    for (var cnt = 0; cnt < coeffs.length; cnt++) {
      out = this.runStage(coeffs[cnt], out);
    }
    return out;
  }

  biquadResponse(FrFsParams params, TempF s) {
    var Fs = params.Fs;
    var Fr = params.Fr;
    // z = exp(j*omega*pi) = cos(omega*pi) + j*sin(omega*pi)
    // z^-1 = exp(-j*omega*pi)
    // omega is between 0 and 1. 1 is the Nyquist frequency.
    var theta = -pi * (Fr / Fs) * 2;
    var z = Complex(cos(theta), sin(theta));
    // k * (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^⁻1 + a2*z^-2)
    var p = s.k * (s.b0 + (z * (s.b1 + s.b2 * z)));
    var q = this.cone + (z * (s.a1 + (s.a2 * z)));
    var h = p / q;
    var res = MagnitudePhase(magnitude: h.magnitude(), phase: h.phase());
    return res;
  }

  calcResponse(FrFsParams params) {
    var res = FilterResponse(
      magnitude: 1,
      phase: 0,
      dBmagnitude: 0,
    );
    for (var cnt = 0; cnt < this.cf.length; cnt++) {
      var r = this.biquadResponse(params, this.cf[cnt]);
      // a cascade of biquads results in the multiplication of H(z)
      // H_casc(z) = H_0(z) * H_1(z) * ... * H_n(z)
      res.magnitude = res.magnitude! * r.magnitude!;
      // phase is wrapped -> unwrap before using
      res.phase = res.phase! + r.phase!;
    }
    res.dBmagnitude = 20 * log(res.magnitude!) * log10e;
    return res;
  }

  reinit() {
    var tempF = List.generate(
        this.f.length,
        (cnt) => TempF(
              b0: Complex(this.f[cnt].b[0], 0),
              b1: Complex(this.f[cnt].b[1], 0),
              b2: Complex(this.f[cnt].b[2], 0),
              a1: Complex(this.f[cnt].a[0], 0),
              a2: Complex(this.f[cnt].a[1], 0),
              k: Complex(this.f[cnt].k, 0),
              z: [0, 0],
            ));
    return tempF;
  }

  calcInputResponse(List<double> input) {
    var tempF = this.reinit();
    return runMultiFilter<List<TempF>>(
        input, tempF, (input, coeffs) => this.doStep(input, coeffs));
  }

  predefinedResponse(double Function(int) def, int length) {
    List<double> input = [];
    for (var cnt = 0; cnt < length; cnt++) {
      input.add(def(cnt));
    }
    final resp = this.calcInputResponse(input);
    var maxFound = false;
    var minFound = false;
    SampleValue? maxSample;
    SampleValue? minSample;
    for (var cnt = 0; cnt < length - 1; cnt++) {
      if (resp[cnt] > resp[cnt + 1] && !maxFound) {
        maxFound = true;
        maxSample = SampleValue(
          sample: cnt,
          value: resp[cnt],
        );
      }

      if (maxFound && !minFound && resp[cnt] < resp[cnt + 1]) {
        minFound = true;
        minSample = SampleValue(
          sample: cnt,
          value: resp[cnt],
        );
        break;
      }
    }
    var ret = Ret(out: resp, min: minSample!, max: maxSample!);
    return ret;
  }

  getComplRes(dynamic n1, dynamic n2) {
    var innerSqrt = pow(n1 / 2, 2) - n2;
    if (innerSqrt < 0) {
      return [
        Complex(-n1 / 2, sqrt(innerSqrt.abs())),
        Complex(-n1 / 2, -sqrt(innerSqrt.abs()))
      ];
    } else {
      return [
        Complex(-n1 / 2 + sqrt(innerSqrt), 0),
        Complex(-n1 / 2 - sqrt(innerSqrt), 0)
      ];
    }
  }

  getPZ() {
    var res = [];
    for (var cnt = 0; cnt < this.cc.length; cnt++) {
      res[cnt] = {};
      res[cnt].z = this.getComplRes(this.cc[cnt].b1, this.cc[cnt].b2);
      res[cnt].p = this.getComplRes(this.cc[cnt].a1, this.cc[cnt].a2);
    }
    return res;
  }

  double singleStep(double input) {
    return this.doStep(input, this.cf);
  }

  multiStep(List<double> input, {bool overwrite = false}) {
    return runMultiFilter<List<TempF>>(
        input, this.cf, (input, coeffs) => this.doStep(input, coeffs),
        overwrite: overwrite);
  }

  filtfilt(List<double> input, {bool overwrite = false}) {
    return runMultiFilterReverse<List<TempF>>(
        runMultiFilter<List<TempF>>(
            input, this.cf, (input, coeffs) => this.doStep(input, coeffs),
            overwrite: overwrite),
        this.cf,
        this.doStep,
        overwrite: true);
  }

  simulate(List<double> input) {
    return this.calcInputResponse(input);
  }

  stepResponse(int length) {
    return this.predefinedResponse((_) => 1.0, length);
  }

  impulseResponse(int length) {
    return this.predefinedResponse((val) {
      if (val == 0) {
        return 1.0;
      } else {
        return 0.0;
      }
    }, length);
  }

  responsePoint(FrFsParams params) {
    return this.calcResponse(params);
  }

  response({int resolution = 100}) {
    var r = resolution * 2;
    var res = List.generate(
        resolution,
        (cnt) => this.calcResponse(FrFsParams(
              Fs: r.toDouble(),
              Fr: cnt.toDouble(),
            )));
    evaluatePhase(res);
    return res;
  }

  polesZeros() {
    return this.getPZ();
  }

  reInit() {
    for (var cnt = 0; cnt < this.cf.length; cnt++) {
      this.cf[cnt].z = [0, 0];
    }
  }
}

class Ret {
  List<double> out;
  SampleValue max;
  SampleValue min;

  Ret({required this.out, required this.min, required this.max});
}

class SampleValue {
  int sample;
  double value;

  SampleValue({required this.sample, required this.value});
}

class FrFsParams {
  double Fr;
  double Fs;

  FrFsParams({required this.Fr, required this.Fs});
}

class TempF {
  Complex b0;
  Complex b1;
  Complex b2;
  Complex a1;
  Complex a2;
  Complex k;
  List<double> z;

  TempF({
    required this.b0,
    required this.b1,
    required this.b2,
    required this.a1,
    required this.a2,
    required this.k,
    required this.z,
  });
}

class CC {
  double a1;
  double a2;
  double b1;
  double b2;

  CC({
    required this.a1,
    required this.a2,
    required this.b1,
    required this.b2,
  });
}

class MagnitudePhase {
  double magnitude;
  double phase;

  MagnitudePhase({required this.magnitude, required this.phase});
}
