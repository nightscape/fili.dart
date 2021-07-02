import 'dart:math';

class FirCoeffs {
  List<double> lowpass(FirParams params) {
    return FirCoeffs.calcImpulseResponse(params);
  }

  List<double> bandpass(FirParams params) {
    return FirCoeffs.invert(FirCoeffs.bs(params));
  }

  List<double> highpass(FirParams params) {
    return FirCoeffs.invert(FirCoeffs.calcImpulseResponse(params));
  }

  List<double> bandstop(FirParams params) {
    return FirCoeffs.bs(params);
  }

  List<double> kbFilter(FirParams params) {
    return FirCoeffs.calcKImpulseResponse(params);
  }

  List<String> available() {
    return ['lowpass', 'highpass', 'bandstop', 'bandpass', 'kbFilter'];
  }

  static double ino(double val) {
    int d = 0;
    double ds = 1;
    double s = 1;
    while (ds > s * 1e-6) {
      d += 2;
      ds *= val * val / (d * d);
      s += ds;
    }
    return s;
  }

  // Kaiser windowd filters
  // desired attenuation can be defined
  // better than windowd sinc filters
  static List<double> calcKImpulseResponse(FirParams params) {
    var Fs = params.Fs!;
    var Fa = params.Fa!;
    var Fb = params.Fb!;
    var o = params.order; //?? 51;
    var alpha = params.Att ?? 100;

    if (o / 2 - (o / 2).floor() == 0) {
      o++;
    }
    var Np = (o - 1) ~/ 2;
    List<double> A = List.filled(Np + 1, 0.0);
    double beta = 0;

    A[0] = 2 * (Fb - Fa) / Fs;
    for (var cnt = 1; cnt <= Np; cnt++) {
      A[cnt] = (sin(2 * cnt * pi * Fb / Fs) - sin(2 * cnt * pi * Fa / Fs)) /
          (cnt * pi);
    }
    // empirical coefficients
    if (alpha < 21) {
      beta = 0;
    } else if (alpha > 50) {
      beta = 0.1102 * (alpha - 8.7);
    } else {
      beta = 0.5842 * pow((alpha - 21), 0.4) + 0.07886 * (alpha - 21);
    }

    var inoBeta = ino(beta);
    List<double> ret = List.filled(Np * 2 + 1, 0.0);
    for (var cnt = 0; cnt <= Np; cnt++) {
      ret[Np + cnt] =
          A[cnt] * ino(beta * sqrt(1 - (cnt * cnt / (Np * Np)))) / inoBeta;
    }
    for (var cnt = 0; cnt < Np; cnt++) {
      ret[cnt] = ret[o - 1 - cnt];
    }
    return ret;
  }

  // note: coefficients are equal to impulse response
  // windowd sinc filter
  static calcImpulseResponse(FirParams params) {
    var Fs = params.Fs!;
    var Fc = params.Fc!;
    var o = params.order;
    var omega = 2 * pi * Fc / Fs;
    double dc = 0;
    List<double> ret = List.filled(o + 1, 0.0);
    // sinc function is considered to be
    // the ideal impulse response
    // do an idft and use Hamming window afterwards
    for (var cnt = 0; cnt <= o; cnt++) {
      if (cnt - o / 2 == 0) {
        ret[cnt] = omega;
      } else {
        ret[cnt] = sin(omega * (cnt - o / 2)) / (cnt - o / 2);
        // Hamming window
        ret[cnt] *= (0.54 - 0.46 * cos(2 * pi * cnt / o));
      }
      dc = dc + ret[cnt];
    }
    // normalize
    for (var cnt = 0; cnt <= o; cnt++) {
      ret[cnt] /= dc;
    }
    return ret;
  }

  // invert for highpass from lowpass
  static List<double> invert(List<double> h) {
    for (var cnt = 0; cnt < h.length; cnt++) {
      h[cnt] = -h[cnt];
    }
    h[(h.length - 1) ~/ 2]++;
    return h;
  }

  static List<double> bs(FirParams params) {
    var lp = FirCoeffs.calcImpulseResponse(FirParams(
      order: params.order,
      Fs: params.Fs,
      Fc: params.F2,
    ));
    var hp = FirCoeffs.invert(FirCoeffs.calcImpulseResponse(FirParams(
      order: params.order,
      Fs: params.Fs,
      Fc: params.F1,
    )));
    List<double> out = [];
    for (var i = 0; i < lp.length; i++) {
      out.add(lp[i] + hp[i]);
    }
    return out;
  }
}

class FirParams {
  int order;
  double? Fa;
  double? Fb;
  double? Fc;
  double? Fs;
  double? F1;
  double? F2;
  double? Att;
  FirParams(
      {required this.order,
      this.Fa,
      this.Fb,
      this.Fc,
      this.Fs,
      this.F1,
      this.F2,
      this.Att});
}
