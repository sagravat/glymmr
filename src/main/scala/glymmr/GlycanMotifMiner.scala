package glymmr

import java.util.Collections

import parser._
import parser.node._
import util.{FalseDiscoveryRate, OneSampleTTest, GlycanInverter}

import scala.Predef
import scala.collection.mutable
import scala.collection.mutable.Map
import scala.collection.parallel.immutable
import scala.collection.parallel.immutable.ParMap
import scala.collection.parallel.mutable.ParSeq
import scala.io.Source



object GlycanMotifMiner {

  val ROOT: Int = 0
  val INTERNAL: Int = 2
  val LEAF: Int = 3


  def getTerminalDeterminants(subsets: List[GlycanNode], glycan: Glycan): List[GlycanNode] =  {

    var motifList: List[GlycanNode] = List[GlycanNode]()

    for (subtree <- subsets) {

      val glycanNode: GlycanNode = glycan.glycanNode
      val foundNode: GNode = glycanNode.find(subtree).get
      val foundNodeChildNodes = foundNode.childList
      val subNodes   = subtree.childList
      if (subtree.toSlc() == foundNode.toSlc()) {
        motifList = subtree :: motifList
      }

      else {
        val nextNode: Option[GNode] = foundNodeChildNodes.find {
          node => !subNodes.exists(x => x.id == node.id) && node.children.size > 0
        }

        if (nextNode.isEmpty) {
          if (glycanNode.find(subNodes.last).get.children.size == 0) {
            motifList = subtree :: motifList

          }
        }
      }
    }

    motifList
  }


  def getInternalMotifMap(subsets: ParMap[GlycanNode, List[GlycanNode]]): Set[String] =  {

    val motifMap: scala.collection.mutable.Map[String,Set[String]] =
      scala.collection.mutable.Map[String,Set[String]]()

    var glycanSlcSet: Set[String] = Set[String]()

    for ( (glycanNode, motifSubsets) <- subsets) {

      for (subtree <- motifSubsets) {

        val foundNode: GNode = glycanNode.find(subtree).get
        val foundNodeChildNodes = foundNode.childList
        val subNodes   = subtree.childList

        if (subtree.toSlc() != foundNode.toSlc())  {

          val nextNode: Option[GNode] = foundNodeChildNodes.find {
            node => !subNodes.exists(x => x.id == node.id) && node.children.size > 0
          }

          if (nextNode.isDefined) {
            val motif: String = subtree.toSlc()
            //            val glycanSlcSet: Set[String] = motifMap.getOrElse(motif, Set())
            //            motifMap += (GlycanInverter.invert(motif) -> (glycanSlcSet + glycanNode.toSlc()) )
            glycanSlcSet += GlycanInverter.invert(motif)
          }

        }
      }
    }
    glycanSlcSet
  }


  def getTerminalMotifs(binders: List[GlycanArrayData],
                        nonBinders: List[GlycanArrayData],
                        minSubsetSize: Int,
                        maxSubsetSize: Int):  List[MotifResult] = {

    // glycanarraydata -> glycan => glycan subtrees

    val nonBinderMotifMap: Map[String, List[GlycanArrayData]] = Map[String, List[GlycanArrayData]]()

    nonBinders.foreach {
      b => {
        // parse code to glycan
        val g: Glycan = GlycanIupacParser.parseToTree(b.linearCode, b.name)

        // get the subsets of the glycan
        val subsets: List[GlycanNode] =
          (for (subsetSize <- minSubsetSize to maxSubsetSize)
            yield compute_subsets(g, subsetSize) ).flatten.toList

        // find the terminal subsets
        val terminalDeterminants: List[GlycanNode] = getTerminalDeterminants(subsets, g)

        // add the motif and list of its glycans to the map
        terminalDeterminants.foreach {

          t => {
            val motif: String = t.toSlc()
            val inverted: String = invert(motif)
            val glycanList: List[GlycanArrayData] = nonBinderMotifMap.getOrElse(inverted, List())
            nonBinderMotifMap += (inverted -> ( b :: glycanList))

          }
        }
      }
    }

    val binderMotifMaps: Map[String, List[GlycanArrayData]] = Map[String, List[GlycanArrayData]]()


    binders.foreach {
      b => {
        // parse code to glycan
        val g: Glycan = GlycanIupacParser.parseToTree(b.linearCode, b.name)

        // get the subsets of the glycan
        val subsets: List[GlycanNode] =
          (for (subsetSize <- minSubsetSize to maxSubsetSize)
            yield compute_subsets(g, subsetSize) ).flatten.toList

        // find the terminal subsets
        val terminalDeterminants: List[GlycanNode] = getTerminalDeterminants(subsets, g)

        // add the motif and list of its glycans to the map
        terminalDeterminants.foreach {

          t => {
            val motif: String = t.toSlc()
            val inverted: String = invert(motif)
            val glycanList: List[GlycanArrayData] = binderMotifMaps.getOrElse(inverted, List())
            binderMotifMaps += (inverted -> ( b :: glycanList))
          }
        }
      }
    }
    val avgBinderRfu: Map[String, Double] = Map[String, Double]()
    var avgBinderRfuList: List[Double] = List[Double]()

    for (
      (motif, refHighList) <- binderMotifMaps;
      refLowList = nonBinderMotifMap.get(motif).getOrElse(List[GlycanArrayData]())
    //        if refHighList.size > refLowList.size
    ) {
      val sum: Double = refHighList.map(r => r.rfuValue).sum.toDouble
      val lowSum: Double = refLowList.map( r => r.rfuValue).sum.toDouble

      val avgVal: Double = (sum+lowSum)/(refHighList.size.toDouble + refLowList.size.toDouble)
      avgBinderRfu += (motif ->  avgVal )
      avgBinderRfuList = avgVal :: avgBinderRfuList
    }

    val binderList: List[(String, List[GlycanArrayData])] = binderMotifMaps.toList

    val pvalues: List[Double] =
      (for ((motif, refHighList) <- binderList;
            pvalue = OneSampleTTest.getLeftTailedPValue(motif, avgBinderRfuList, avgBinderRfu.getOrElse(motif, 0.0))
      ) yield pvalue)

    val fdrResult: FDR_Result = FalseDiscoveryRate.fdr(pvalues, 0.05)

    val motifResults: List[MotifResult] =
      for (((motif, refHighList),index) <- binderList.zipWithIndex;
           refLowList = nonBinderMotifMap.get(motif).getOrElse(List[GlycanArrayData]());
           pvalue = OneSampleTTest.getLeftTailedPValue(motif, avgBinderRfuList, avgBinderRfu.getOrElse(motif, 0.0));
           qvalue = fdrResult.adj_p(index)
           //        if refHighList.size > refLowList.size && pvalue <= 0.05
           if pvalue <= 0.05
      )
        yield
          MotifResult(
            motif,
            refHighList,
            refLowList,
            avgBinderRfu.getOrElse(motif, 0.0),
            pvalue,
            fdrResult.crit_p,
            qvalue
          )

    //      motifResults.filter( m => m.refHigh.size > m.refLow.size).toList.sortBy( m => m.motif.size)
    motifResults.toList.sortBy(m => m.pValue )

  }


  def main(args: Array[String]) {

    val file:String = args(0)
    val subsetSize: Int = args(1).toInt
    val motifType: String = args(2)
    var start: Long = 0l
    var stop: Long = 0l
    val trees: Seq[Glycan] =
      Source.fromFile(file).getLines.
        toIndexedSeq.map(line => GlycanIupacParser.parseToTree(line + "-Sp0"))

    (2 to 3).foreach(
      subsetSize => {

        start = System.currentTimeMillis / 1000
        // subsets
        val subsets: Predef.Map[Glycan, List[GlycanNode]] =
          trees.map( glycan => glycan -> compute_subsets(glycan, subsetSize)).toMap
        stop = System.currentTimeMillis / 1000

        if (motifType == "terminal") {
          //          var motifs: Set[String] = Set[String]()
          var motifs: List[String] = List[String]()
          subsets.keys.foreach( (glycanKey) => {
            val terminalMotifs: List[GlycanNode] =
              getTerminalDeterminants(subsets.get(glycanKey).get, glycanKey)
            terminalMotifs.foreach( (motif) => {
              if (!motifs.contains(motif.toSlc())) {

                motifs = motif.toSlc() :: motifs
              }
            })

          }
          )
          println(s"${file}, ${subsetSize}, ${motifs.size}")
          for ( motif <- motifs )
          {
            println(GlycanInverter.invert(motif))
          }


          //          val terminalMotifs: Map[String,List[String]] = getTerminalMotifs(subsets, List(), 2, 4)
          //
          //          for ( (motif, glycanSet) <- terminalMotifs) {
          //            println(s"${motif},${glycanSet.size}")
          //          }
        }
        else {

          //          val internalMotifMap: Set[String] = getInternalMotifMap(subsets)
          //
          //          for ( (motif) <- internalMotifMap) {
          //            println(s"${motif}")
          //          }

        }

      }
    )

  }



  def compute_subsets(glycan: Glycan, subsetSize: Int): List[GlycanNode] = {

    val tree: GlycanNode = glycan.glycanNode
    //    log.info(tree.toSlc() + " subset size: " + subsetSize)
    var validSubsets: List[GlycanNode] = List[GlycanNode]()
    val children = tree.childList
    for (nodes <- children.combinations(subsetSize)) {
      var isCorrect = true
      var i = 0
      var subtree: GlycanNode = null
      //println(subsetSize + ", " + nodes.size + ", " + nodes)
      for (glycanNode <- nodes) {

        glycanNode match {

          case g: GlycanNode => {
            val parentId = glycanNode.parentId
            if (parentId != 0 && i > 0 && isCorrect) {
              if (nodes.filter(_.id == parentId).size == 0) {
                isCorrect = false
              }
              else {
                val node =
                  GlycanNode(g.id, g.slc,
                    g.anomer,
                    g.link,
                    List())(parentId, if (g.children.size == 0) LEAF else INTERNAL)
                subtree.findParent(node) match {
                  case Some(GlycanNode(_, _, _, _, _)) => {
                    subtree.findParent(node).get.add(node)

                  }
                  case None => {
                    //println("stop!")
                  }
                }
              }
            }
            else {
              subtree = new GlycanNode(g.id,
                g.slc,
                g.anomer,
                g.link,
                List())(parentId, g.nodeType)
            }
            i = i + 1

          }
        }
      }
      if (isCorrect) validSubsets = subtree :: validSubsets
    }


    validSubsets
  }

  def getFlattenedMotifsFromList(subsets: List[GlycanNode],
                                 nodeTypeFilter: GNode => Boolean): List[GlycanNode] = {

    subsets.filter(_.childList.exists(nodeTypeFilter)).distinct

  }

  def getFlattenedMotifs(subsetMap: Map[String, List[GlycanNode]],
                         nodeTypeFilter: GNode => Boolean): Iterable[List[GlycanNode]] = {

    for {
      (key,listNodes) <- subsetMap
    }
      yield listNodes.filter(_.childList.exists(nodeTypeFilter)).distinct

    //subsetMap.values.filter(_.exists(nodeTypeFilter)).toList.flatten.distinct
  }



  def invert(s: String): String = {

    if (s.length < 3) return s
    val residues: java.util.List[String] = new java.util.ArrayList[String]

    var first: String = ""
    var rest: String = s


    if ((s.charAt(2)+"").matches("[A-Z]") || s.charAt(2) == '(' ) {
      first = s.substring(0,2)
      rest = s.substring(2)
    }


    val length: Int = rest.length
    var isParen: Boolean = (rest.startsWith("(") || first.startsWith("("))
    //for (int i = 0; i < length;) {
    while (!(rest == "")) {
      var residue: String = ""
      if (isParen) {
        val paren: String = rest.substring(0, 1)
        residues.add(if ((paren == ")")) "(" else ")")
      }
      else {
        try {

          residue = rest.substring(0, 3)
        } catch {
          case e: Exception =>
            println(s)
        }
        residues.add(residue)
      }
      if (isParen) {
        rest = rest.substring(1)
        isParen = false
      }
      else {
        rest = rest.substring(3)
      }
      if (rest.startsWith("(") || rest.startsWith(")")) {
        isParen = true
      }
      else {
      }
    }
    Collections.reverse(residues)
    val sb: StringBuilder = new StringBuilder
    import scala.collection.JavaConversions._
    for (residue <- residues) {
      sb.append(residue)
    }
    sb.append(first)

    return sb.toString
  }

}
