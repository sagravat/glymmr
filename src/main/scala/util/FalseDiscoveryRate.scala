package util

import glymmr.FDR_Result

/**
  * Created by sagravat on 2/21/16.
  */
object FalseDiscoveryRate {


  def fdr(pvalues: List[Double], q: Double): FDR_Result = {

    val p_sorted: List[Double] = pvalues.sorted
    val sort_ids : Array[Int] = pvalues.zipWithIndex.sortBy(_._1).map(_._2).toArray
    val unsorted_ids : Array[Int] = sort_ids.zipWithIndex.sortBy(_._1).map(_._2)

    var m: Int = pvalues.size
    val thresh: List[Double] = (1 to m).map( i => (i * q)/m).toList

    val wtd_p: List[Double] = p_sorted.map( p => p * m ).zipWithIndex.map( e => e._1/(e._2+1) )

    val wtd_p_sorted: List[Double] = wtd_p.sorted
    val wtd_p_sindex: Array[Int] = wtd_p.zipWithIndex.sortBy(_._1).map(_._2).toArray
    val adj_p: Array[Double] = Array.fill[Double](m)(0.0)
    var nextfill = 1
    val length = m
    for (k <- 0 until m) {
      if (wtd_p_sindex(k) >= nextfill) {
        for (j <- nextfill to wtd_p_sindex(k)) {
          adj_p(j) = wtd_p_sorted(k)
        }
        nextfill = wtd_p_sindex(k) + 1
        if (nextfill > m) {
          m = Int.MaxValue
        }
      }
    }

    val new_adj: Array[Double] = Array.fill[Double](m)(0.0)
    for (k <- 0 until length) {

      new_adj(k) = adj_p( unsorted_ids(k) )

    }

    val zipped: List[(Double, Double)] = p_sorted.zip(thresh)
    val rejected: List[(Double, Double)] = zipped.filter( f => f._1 <= f._2 )
    val max_id: Int = rejected.size - 1



    if (max_id <= 0) {

      FDR_Result(Array.fill[Int](m)(0).toList, 0.0, new_adj.toList)

    }
    else {
      val crit_p: Double = p_sorted(max_id)
      val h: List[Int] = pvalues.map( p => if (p <= crit_p) 1 else 0)

      FDR_Result(h, crit_p, new_adj.toList)
    }

  }

}
