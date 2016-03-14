package util

import org.apache.commons.math3.special.Beta
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.stat.inference.MannWhitneyUTest
import org.apache.commons.math3.util.FastMath

/**
 * Created by Sanjay Agravat on 8/28/15.
 */
object OneSampleTTest {

  def main(args: Array[String]) {



  }

  def mannWhitneyUTest(motif: String, x: List[Double], y: List[Double]): Double = {


    val xArr: Array[Double] = x.toArray
    val yArr: Array[Double] = y.toArray

    val pvalue: Double = new MannWhitneyUTest().mannWhitneyUTest(xArr, yArr)
    println(motif + "\n\t" + x + "\n\t" + y + "\n\t" + pvalue)

    pvalue

  }

  def getLeftTailedPValue(motif: String, rfuValues: List[Double], mu: Double): Double = {
    
    val samples: Array[Double] = rfuValues.toArray

    val tStatistic: Double =
      (t(StatUtils.mean(samples),
        mu,
        StatUtils.variance(samples),
        samples.length)
      )
    val pvalue: Double = cumulativeProbability(tStatistic, samples.length - 1)

    println(s"${motif} mu ${mu}, t=${tStatistic}, p=${pvalue}")

    pvalue
  }


  def t(m: Double, mu: Double, v: Double, n: Double): Double = {
    (m - mu) / FastMath.sqrt(v / n)
  }


  def cumulativeProbability(x: Double, degreesOfFreedom: Int): Double = {
    var ret: Double = 0.0
    if (x == 0) {
      ret = 0.5
    } else {
      val t:Double  =
        Beta.regularizedBeta(
          degreesOfFreedom / (degreesOfFreedom + (x * x)),
          0.5 * degreesOfFreedom,
          0.5);
      if (x < 0.0) {
        ret = 0.5 * t;
      }
      else if (x < 1.0) {
        ret = t;
      }
      else {
        ret = 1.0 - 0.5 * t;

      }
    }

    return ret;
  }

}
