import 'dart:math';

import 'package:fili.dart/fft.dart';
import 'package:test/test.dart';

main() {
  List<double> sinewave = [];

  final fftCalc = Fft(8192);
  for (var cnt = 0; cnt < 8192; cnt++) {
    sinewave.add(sin(2 * pi * (113.33232 * cnt / 8192)));
  }

  group('fft_calc', () {
    final fftResult = fftCalc.forward(sinewave, 'hanning');
    final magnitude = fftCalc.magnitude(fftResult);
    test('can forward', () {
      expect(fftResult.re.length, equals(8192));
      expect(fftResult.im.length, equals(8192));
    });

    test('can calculate the magnitude', () {
      expect(magnitude.length, equals(8192));
    });

    test('can calculate the phase', () {
      var phase = fftCalc.phase(fftResult);
      expect(phase.length, equals(8192));
    });

    test('can calculate the magnitude in dB', () {
      var magdb = fftCalc.magToDb(magnitude);
      expect(magdb.length, equals(8192));
    });

    test('can inverse', () {
      var original = fftCalc.inverse(fftResult);
      expect(original.length, equals(8192));
    });

    test('can calculate with all window functions', () {
      var w = fftCalc.windows();
      for (var i = 0; i < w.length; i++) {
        final fftResult = fftCalc.forward(sinewave, w[i]);
        expect(fftResult.re.length, equals(8192));
        expect(fftResult.im.length, equals(8192));
      }
    });
  });

  group('fft_helpers', () {
    test('can get available windows', () {
      var w = fftCalc.windows();
      expect(w.length, isNot(equals(0)));
    });

    test('detects wrong radix', () {
      expect(
          () => Fft(1234),
          throwsA(const TypeMatcher<Exception>()
              .having((e) => e.toString(), 'text', matches(".*power of 2.*"))));
    });
  });
}
