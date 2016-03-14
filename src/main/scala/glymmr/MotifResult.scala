package glymmr

/**
 * Created by Sanjay Agravat on 8/22/15.
 */
case class MotifResult(motif: String,
                       refHigh: List[GlycanArrayData],
                       refLow: List[GlycanArrayData],
                       avgRfu: Double=0.0,
                       pValue: Double=0.0,
                       criticalValue: Double=0.0,
                       qValue: Double=0.0) {

}
