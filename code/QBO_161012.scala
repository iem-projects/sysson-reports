// version: 12-Oct-2016

implicit class MyGEOps(private val in: GE) /* extends AnyVal */ {
  def +: (head: Constant): GE = Flatten(Seq(head, in))
  def :+ (last: Constant): GE = Flatten(Seq(in, last))
  def ++ (that: GE      ): GE = Flatten(Seq(in, that))
  def \  (r: Range      ): GE = r.map(i => in \ i): GE
}

// ---- matrix ----

val vr    = Var("anom")
val dTime = Dim(vr, "time")
val dAlt  = Dim(vr, "altitude")

val speed = UserValue.kr("speed", 6)
val tp    = dTime.play(speed)
val timeTr= Impulse.ar(speed)
val vp    = vr.play(tp, interp = 1 /* 2 */)

// ---- configuration ----

// maximum number of trajectories followed
val numTraj       = 4

// maximum number of concurrently playing voices
val numVoices     = numTraj * 2

// maximum jump in altitude (normalized - full range equals one) before trajectory is interrupted
val maxDf         = UserValue.kr("traj-max-dif [0-1]", 0.25)

// trajectory amplitude envelope attack duration
val egAtk         = UserValue.kr("traj-atk [s]", 0.2)

// trajectory amplitude envelope release duration
val egRls         = UserValue.kr("traj-rls [s]", 0.5)

// trajectory frequency and amplitude smear time
val lagTime       = 1.0

// threshold above which temperatures are considered anormal
val magThresh     = UserValue.kr("mag-thresh [°]", 1.5)

// maximum magnitude considered for anormal values
val magMax        = UserValue.kr("mag-max [°]", 8.0)

// sonification minimum oscillator frequency for hot anomalies
val minFreqH      = UserValue.kr("min-freq hot [Hz]", 300)

// sonification maximum oscillator frequency for hot anomalies
val maxFreqH      = UserValue.kr("max-freq hot [Hz]", 600)

// sonification minimum oscillator frequency for cold anomalies
val minFreqC      = UserValue.kr("min-freq cold [Hz]", 1000)

// sonification maximum oscillator frequency for cold anomalies
val maxFreqC      = UserValue.kr("max-freq cold [Hz]", 2000)

// sonification maximum oscillator frequency for cold anomalies
val filterQ       = UserValue.kr("filter Q (1 to 100)", 20)
val filterRQ      = filterQ.reciprocal

// time grid indicator volume
val gridAmp       = UserValue.kr("time grid [dB]", -18).dbamp

// ---- state ----

val stateInKr     = LocalIn.kr(Seq.fill(numVoices * 3 * 2)(0))

// ---- time grid ----

def mkGrid(): GE = {
  // XXX TODO --- we need an operator
  // that can automatically decode the dates,
  // since other data sets will use other scaling
  val seconds      = tp
  val daysPerYear  = 365.2422  // average according to NASA
  val daysPerMonth = daysPerYear / 12
  val secsPerMonth = daysPerMonth * 24 * 60 * 60
  val months       = seconds / secsPerMonth
  val month        = months.floor % 12 // + 1
  val isJan        = month sig_== 0
  val isJul        = month sig_== 6
  
  val monthPulse  = timeTr
  val gridDecTime = speed.reciprocal.min(0.5)
  val gridDecay   = Decay.ar(monthPulse, gridDecTime)
  val gridBase    = WhiteNoise.ar(gridDecay) * gridAmp
  val gridJan     = Resonz.ar(gridBase, 1000, rq = 1)
  val gridJul     = Resonz.ar(gridBase, 3500, rq = 0.5) * 2 // * LFPulse.ar(40)
  val gridPlain0  = Resonz.ar(gridBase, 2000, rq = 2) * 0.333
  val gridPlain   = Select.ar(speed <= 6, Seq(DC.ar(0), gridPlain0))
  val gridIdx     = isJan | (isJul << 1)
  val gridSig     = Select.ar(index = gridIdx, in = Seq(gridPlain, gridJan, gridJul))
  
  Pan2.ar(gridSig)
}

// ---- matrix analysis ----

val voiceNos      = 0 until numVoices: GE
val numAlt        = dAlt.size
val maskWidth     = numAlt / 2 // 4
val vpChanIdx     = ChannelIndices(vp)

def mkSide(isUp: Boolean): (GE, GE) = {
  // ---- state ----
  val stateOff      = if (isUp) 0 else 3
  var voiceFreq     = stateInKr \ ((numVoices * (stateOff+0)) until (numVoices * (stateOff+1)))
  var voiceAmp      = stateInKr \ ((numVoices * (stateOff+1)) until (numVoices * (stateOff+2)))
  var voiceOnOff    = stateInKr \ ((numVoices * (stateOff+2)) until (numVoices * (stateOff+3)))
  
  // ---- trace trajectories ----
  
  def extract(in: GE, res: Seq[(GE, GE)], trjIdx: Int): Seq[(GE, GE)] = 
    if (trjIdx == numTraj) res else {
      val (bestIdx, bestVal) = 
        if (isUp) {
          val best = ArrayMax.kr(in)
          best.index -> best.value
        } else {
          val best = ArrayMin.kr(in)
          best.index -> best.value
        }

      val freq0     = bestIdx / (numAlt - 1) // .linexp(0, numAlt - 1, 200, 4000)
      val amp0      = if (isUp)
        bestVal.clip(magThresh, magMax).linlin(magThresh, magMax, 0, 1)
      else
        bestVal.clip(-magMax, -magThresh).linlin(-magMax, -magThresh, 1, 0)
      
      val mask      = in * vpChanIdx.absdif(bestIdx) > maskWidth
      extract(in = mask, res = res :+ (freq0 -> amp0), trjIdx = trjIdx + 1)
    }
  
  var activated   = Vector.fill(numVoices)(0: GE): GE
  
  val (freqInSq, ampInSq) = extract(in = A2K.kr(vp), res = Vector.empty, trjIdx = 0).unzip
  val freqIn = freqInSq: GE
  val ampIn  = ampInSq : GE
  
  // for each frequency, find the best past match
  val noFounds = (0 until numTraj).map { tIdx =>
    val fIn         = freqIn \ tIdx
    val aIn         = ampIn  \ tIdx
    val isOn        = aIn > 0
  
    val freqMatch   = (maxDf - (voiceFreq absdif fIn)).max(0)
    val bothOn      = voiceOnOff & isOn
    val bestIn      = 0 +: (freqMatch * (bothOn & !activated))
    val best        = ArrayMax.kr(bestIn)
    val bestIdx     = best.index - 1
  
    val bestMask    = voiceNos sig_== bestIdx
    activated      |= bestMask
    val bestMaskN   = !bestMask
    voiceFreq       = voiceFreq * bestMaskN + fIn * bestMask
    voiceAmp        = voiceAmp  * bestMaskN + aIn * bestMask
    
    bestIdx sig_== -1
  }
  
  for (tIdx <- 0 until numTraj) {
    val fIn             = freqIn \ tIdx
    val aIn             = ampIn  \ tIdx
    val isOn            = aIn > 0
    val voiceAvail      = !(activated | voiceOnOff)
  
    val notFound        = noFounds(tIdx)
    val startTraj       = notFound & isOn
    val free            = ArrayMax.kr(0 +: (startTraj & voiceAvail))
    val freeIdx         = free.index - 1
    val freeMask        = voiceNos sig_== freeIdx
    activated          |= freeMask
    val freeMaskN       = !freeMask
    voiceFreq           = voiceFreq * freeMaskN + fIn * freeMask
    voiceAmp            = voiceAmp  * freeMaskN + aIn * freeMask
  }
  
  // ---- voice generation ----
  val voiceEnv      = Env.asr(attack = egAtk, release = egRls)
  val voiceEG       = EnvGen.ar(voiceEnv, gate = activated)
  
  // ---- state out ----
  val voiceEGOn = A2K.kr(voiceEG) sig_!= 0
  voiceOnOff    = activated | voiceEGOn
  
  val stateOutKr  = voiceFreq ++ voiceAmp ++ voiceOnOff
  
  // ---- sound generation ----
  
  // gate so attack doesn't lag
  val lagTimeGt = (activated sig_== Delay1.kr(activated)) * lagTime
  val ampScale  = voiceAmp.linexp(0, 1, -20.dbamp, 0.dbamp)
  val ampLag    = Lag.ar(voiceAmp , time = lagTimeGt)
  val minFreq   = if (isUp) minFreqH else minFreqC
  val maxFreq   = if (isUp) maxFreqH else maxFreqC
  val freqScale = voiceFreq.linexp(0, 1, minFreq, maxFreq)
  val freqLag   = Lag.ar(freqScale, time = lagTimeGt)
  val osc       = if (isUp) {
      LPF.ar(Saw.ar(freqLag/2), freqLag) * 0.25 +
      Resonz.ar(BrownNoise.ar.sqrt, freqLag, rq = filterRQ) * 1.5 // * 1
    } else {
      Resonz.ar(Dust2.ar(400), freqLag, rq = filterRQ) * 10
    }
  val sines     = osc * voiceEG * ampLag
  val mix       = Mix(sines)
  (mix -> stateOutKr)
}

val (left , leftState)  = mkSide(true )
val (right, rightState) = mkSide(false)

LocalOut.kr(leftState ++ rightState)

val grid  = mkGrid()
val sig   = (Seq(left, right): GE) + grid
val mix   = Limiter.ar(sig)

output := mix

