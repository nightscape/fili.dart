import 'package:fili.dart/fir_coeffs.dart';
import 'package:fili.dart/fir_filter.dart';
import 'package:test/test.dart';

main() {
  var firCalculator = new FirCoeffs();

  group('fir-lp', () {
    List<double> filterCoeffs;

    filterCoeffs =
        firCalculator.lowpass(FirParams(order: 100, Fs: 4000, Fc: 500));
    FirFilter filter = new FirFilter(filterCoeffs);
    test('can calculate coeffs', () {
      // // filterCoeffs.should.be.an.Array;
      // // filterCoeffs[44].should.be.a.Number;
      expect(filterCoeffs.length, equals(101));
    });

    test('can generate a filter', () {
      // // filter.should.be.an.Object;
    });

    test('can do a single step', () {
      var out = filter.singleStep(10);
      // // out.should.be.a.Number;
      expect(out, isNot(equals(0)));
    });

    test('can do multiple steps', () {
      List<double> simInput = [];
      for (var i = 0; i < 10000; i++) {
        simInput.add(i % 10 - 5);
      }
      ;
      var out = filter.multiStep(simInput);
      // // out.should.be.an.Array;
      expect(out.length, equals(10000));
      expect(out[111], isNot(equals(simInput[111])));
    });

    test('can simulate multiple steps', () {
      List<double> simInput = [];
      for (var i = 0; i < 10000; i++) {
        simInput.add(i % 10 - 5);
      }
      ;
      var out = filter.simulate(simInput);
      // // out.should.be.an.Array;
      expect(out.length, equals(10000));
      expect(out[111], isNot(equals(simInput[111])));
    });

    test('calculates filter response', () {
      var r = filter.response(200);
      // // r.should.be.an.Array;
      expect(r.length, equals(200));
      // // r[20].should.be.an.Object;
      // // r[20].magnitude.should.be.a.Number;
      // // r[20].dBmagnitude.should.be.a.Number;
      // // r[20].phase.should.be.a.Number;
      // // r[20].unwrappedPhase.should.be.a.Number;
      // // r[20].phaseDelay.should.be.a.Number;
      // // r[20].groupDelay.should.be.a.Number;

      expect(r[20].magnitude, isNot(equals(0)));
      expect(r[20].dBmagnitude, isNot(equals(0)));
      expect(r[20].phase, isNot(equals(0)));
      expect(r[20].unwrappedPhase, isNot(equals(0)));
      expect(r[20].phaseDelay, isNot(equals(0)));
      expect(r[20].groupDelay, isNot(equals(0)));
    });

    test('reinit does not crash', () {
      filter.reinit();
    });
  });

  group('fir-hp', () {
    List<double> filterCoeffs =
        firCalculator.highpass(FirParams(order: 100, Fs: 4000, Fc: 1457));
    FirFilter filter = new FirFilter(filterCoeffs);
    test('can calculate coeffs', () {
      expect(filterCoeffs.length, equals(101));
    });

    test('can do a single step', () {
      var out = filter.singleStep(10);
      expect(out, isNot(equals(0)));
    });

    test('can do multiple steps', () {
      List<double> simInput = [];
      for (var i = 0; i < 10000; i++) {
        simInput.add(i % 10 - 5);
      }
      var out = filter.multiStep(simInput);
      expect(out.length, equals(10000));
      expect(out[111], isNot(equals(simInput[111])));
    });

    test('can simulate multiple steps', () {
      List<double> simInput = [];
      for (var i = 0; i < 10000; i++) {
        simInput.add(i % 10 - 5);
      }
      var out = filter.simulate(simInput);
      expect(out.length, equals(10000));
      expect(out[111], isNot(equals(simInput[111])));
    });

    test('calculates filter response', () {
      var r = filter.response(200);
      expect(r.length, equals(200));

      expect(r[20].magnitude, isNot(equals(0)));
      expect(r[20].dBmagnitude, isNot(equals(0)));
      expect(r[20].phase, isNot(equals(0)));
      expect(r[20].unwrappedPhase, isNot(equals(0)));
      expect(r[20].phaseDelay, isNot(equals(0)));
      expect(r[20].groupDelay, isNot(equals(0)));
    });

    test('reinit does not crash', () {
      filter.reinit();
    });
  });

  group('fir-br', () {
    List<double> filterCoeffs = firCalculator
        .bandstop(FirParams(order: 100, Fs: 4000, F1: 457, F2: 1457));
    FirFilter filter = new FirFilter(filterCoeffs);
    test('can calculate coeffs', () {
      expect(filterCoeffs.length, equals(101));
    });

    test('can do a single step', () {
      var out = filter.singleStep(10);
      expect(out, isNot(equals(0)));
    });

    test('can do multiple steps', () {
      List<double> simInput = [];
      for (var i = 0; i < 10000; i++) {
        simInput.add(i % 10 - 5);
      }
      ;
      var out = filter.multiStep(simInput);
      expect(out.length, equals(10000));
      expect(out[111], isNot(equals(simInput[111])));
    });

    test('can simulate multiple steps', () {
      List<double> simInput = [];
      for (var i = 0; i < 10000; i++) {
        simInput.add(i % 10 - 5);
      }
      ;
      var out = filter.simulate(simInput);
      expect(out.length, equals(10000));
      expect(out[111], isNot(equals(simInput[111])));
    });

    test('calculates filter response', () {
      var r = filter.response(200);
      expect(r.length, equals(200));

      expect(r[20].magnitude, isNot(equals(0)));
      expect(r[20].dBmagnitude, isNot(equals(0)));
      expect(r[20].phase, isNot(equals(0)));
      expect(r[20].unwrappedPhase, isNot(equals(0)));
      expect(r[20].phaseDelay, isNot(equals(0)));
      expect(r[20].groupDelay, isNot(equals(0)));
    });

    test('reinit does not crash', () {
      filter.reinit();
    });
  });

  group('fir-bp', () {
    List<double> filterCoeffs = firCalculator
        .bandpass(FirParams(order: 100, Fs: 4000, F1: 577, F2: 1111));
    FirFilter filter = new FirFilter(filterCoeffs);

    test('can calculate coeffs', () {
      expect(filterCoeffs.length, equals(101));
    });

    test('can do a single step', () {
      var out = filter.singleStep(10);
      expect(out, isNot(equals(0)));
    });

    test('can do multiple steps', () {
      List<double> simInput = [];
      for (var i = 0; i < 10000; i++) {
        simInput.add(i % 10 - 5);
      }
      ;
      var out = filter.multiStep(simInput);
      expect(out.length, equals(10000));
      expect(out[111], isNot(equals(simInput[111])));
    });

    test('can simulate multiple steps', () {
      List<double> simInput = [];
      for (var i = 0; i < 10000; i++) {
        simInput.add(i % 10 - 5);
      }
      ;
      var out = filter.simulate(simInput);
      expect(out.length, equals(10000));
      expect(out[111], isNot(equals(simInput[111])));
    });

    test('calculates filter response', () {
      var r = filter.response(200);
      expect(r.length, equals(200));

      expect(r[20].magnitude, isNot(equals(0)));
      expect(r[20].dBmagnitude, isNot(equals(0)));
      expect(r[20].phase, isNot(equals(0)));
      expect(r[20].unwrappedPhase, isNot(equals(0)));
      expect(r[20].phaseDelay, isNot(equals(0)));
      expect(r[20].groupDelay, isNot(equals(0)));
    });

    group('fir-kb-bp', () {
      List<double> filterCoeffs = firCalculator.kbFilter(
          FirParams(order: 101, Fs: 4000, Fa: 577, Fb: 1111, Att: 100));
      FirFilter filter = new FirFilter(filterCoeffs);

      test('can calculate coeffs', () {
        expect(filterCoeffs.length, equals(101));
      });

      test('can do a single step', () {
        var out = filter.singleStep(10);
        expect(out, isNot(equals(0)));
      });

      test('can do multiple steps', () {
        List<double> simInput = [];
        for (var i = 0; i < 10000; i++) {
          simInput.add(i % 10 - 5);
        }
        var out = filter.multiStep(simInput);
        expect(out.length, equals(10000));
        expect(out[111], isNot(equals(simInput[111])));
      });

      test('can simulate multiple steps', () {
        List<double> simInput = [];
        for (var i = 0; i < 10000; i++) {
          simInput.add(i % 10 - 5);
        }
        var out = filter.simulate(simInput);
        expect(out.length, equals(10000));
        expect(out[111], isNot(equals(simInput[111])));
      });

      test('calculates filter response', () {
        var r = filter.response(200);
        expect(r.length, equals(200));

        expect(r[20].magnitude, isNot(equals(0)));
        expect(r[20].dBmagnitude, isNot(equals(0)));
        expect(r[20].phase, isNot(equals(0)));
        expect(r[20].unwrappedPhase, isNot(equals(0)));
        expect(r[20].phaseDelay, isNot(equals(0)));
        expect(r[20].groupDelay, isNot(equals(0)));
      });
    });

    test('reinit does not crash', () {
      filter.reinit();
    });
  });

  group('fir-helpers', () {
    test('can get available filters', () {
      var av = firCalculator.available();
      expect(av.length, isNot(equals(0)));
      // // av[1].should.be.a.String;
    });
  });
}
