;; 1. Author: R.C. van Wijk (r.c.van.wijk@lacdr.leidenuniv.nl)
;; 3. Label: Final metabolite model
;; 4. Dataset: Paracetamol, paracetamol-glucuronide, paracetamol-sulphate in zebrafish larvae of 5 days post fertilization in homogenates, blood samples, and excreted into medium.

$PROBLEM    PK
$INPUT      ID TIME AMT DV EVID MDV CMT XEXP BQL XEXP2 AGE BLOOD
$DATA      Real_Zebrafish_paracetamol_metabolite.csv IGNORE=@ IGNORE=(BQL.GT.0)

; units
; TIME = min
; DV = pmole / larva or pmole / nL
; V = nL
; kA = pmole / min

$SUBROUTINE ADVAN13 TOL=9
$MODEL      COMP ; CMT 1 dosing compartment
            COMP ; CMT 2 central paracetamol in larva
            COMP ; CMT 3 central paracetamol-glucuronide in larva
            COMP ; CMT 4 central paracetamol-sulphate in larva
            COMP ; CMT 5 excreted paracetamol in medium
            COMP ; CMT 6 excreted paracetamol-glucuronide in medium
            COMP ; CMT 7 excreted paracetamol-sulphate in medium
$PK      
TVK25 = THETA(1)
TVK12 = THETA(2)
TVK23 = THETA(3)/1000
TVK24 = THETA(4)
TVK36 = THETA(5)/1000
TVK47 = THETA(6)/1000
TVV2 = THETA(7)
TVV3 = THETA(8)
TVV4 = THETA(9)
TVT50 = THETA(10)
TVF5 = THETA(11)
TVF6 = THETA(12)/1000
TVF7 = THETA(13)

K25 = TVK25 * EXP(ETA(1)) ;first order rate of excretion of paracetamol
K12 = TVK12               ;zero order rate of absorption
K23 = TVK23               ;first order metabolic formation rate for glucuronidation
K24 = TVK24               ;time-dependent metabolic formation rate for sulphation
K36 = TVK36               ;first order rate of excretion of paracetamol-glucuronide
K47 = TVK47               ;first order rate of excretion of paracetamol-sulphate       
V2 = TVV2                 ;distribution volume paracetamol
V3 = TVV3                 ;distribution volume paracetamol-glucuronide
V4 = TVV4                 ;distribution volume paracetamol-sulphate
T50 = TVT50               ;time at which 50% of maximal sulphation is reached
F5 = TVF5                 ;fraction excreted parent retrieved from medium sample
F6 = TVF6                 ;fraction excreted glucuronide-metabolite retrieved from medium sample
F7 = TVF7                 ;fraction excreted sulphate-metabolite retrieved from medium sample 

ET1 = ETA(1)

S2 = V2
S3 = V3
S4 = V4

$DES      
DADT(1) = 0
DADT(2) = K12 * A(1) - K25 * A(2) - K23 * A(2) - K24 * (1 - T/(T50+T)) * A(2)
DADT(3) = K23 * A(2) - K36 * A(3)
DADT(4) = K24 * (1 - T/(T50+T)) * A(2) - K47 * A(4)
DADT(5) = K25 * A(2) * F5
DADT(6) = K36 * A(3) * F6
DADT(7) = K47 * A(4) * F7

$ERROR      
IF(CMT.EQ.2.AND.BLOOD.EQ.0) THEN
IPRED = A(2) 
Y = IPRED * (1 + EPS(1)) + EPS(7) ; comb error paracetamol in larva
ENDIF
IF(CMT.EQ.3.AND.BLOOD.EQ.0) THEN
IPRED = A(3)
Y = IPRED * (1 + EPS(2)) + EPS(8) ; comb error paracetamol-glucuronide in larva
ENDIF
IF(CMT.EQ.4.AND.BLOOD.EQ.0) THEN
IPRED = A(4) 
Y = IPRED * (1 + EPS(3)) + EPS(9) ; prop error paracetamol-sulphate in larva
ENDIF

IF(CMT.EQ.2.AND.BLOOD.EQ.1) THEN
IPRED = A(2)/V2 
Y = IPRED * (1 + EPS(4)) + EPS(10) ; add error paracetamol in blood
ENDIF
IF(CMT.EQ.3.AND.BLOOD.EQ.1) THEN
IPRED = A(3)/V3
Y = IPRED * (1 + EPS(5)) + EPS(11) ; prop error paracetamol-glucuronide in blood
ENDIF
IF(CMT.EQ.4.AND.BLOOD.EQ.1) THEN
IPRED = A(4)/V4 
Y = IPRED * (1 + EPS(6)) + EPS(12) ; prop error paracetamol-sulphate in blood
ENDIF

IF(CMT.EQ.5) THEN
IPRED = A(5)
Y = IPRED * (1 + EPS(13)) + EPS(16) ; add error paracetamol in medium
ENDIF
IF(CMT.EQ.6) THEN
IPRED = A(6)
Y = IPRED * (1 + EPS(14)) + EPS(17) ; prop error paracetamol-glucuronide in medium
ENDIF
IF(CMT.EQ.7) THEN
IPRED = A(7)
Y = IPRED * (1 + EPS(15)) + EPS(18) ; comb error paracetamol-sulphate in medium
ENDIF

IRES = DV - IPRED

$THETA  (0,0.01142) ; K25
 (0,0.3990)         ; K12
 (0,2.592)          ; K23
 (0,0.09168)        ; K24
 (0,5.838)          ; K36
 (0,0.02447)        ; K47
 (0,103)            ; V2
 (0,240)            ; V3
 (0,262)            ; V4
 (0,5.253)          ; T50
 (0, 0.01, 1)       ; f_Q
 (0, 0.01, 1)       ; f_G
 (0, 0.2, 1)        ; f_S
 
$OMEGA  0  FIX      ; IIV K25

$SIGMA  BLOCK(3)
 0.1366             ; prop error paracetamol in larva
 -0.1348 0.3150     ; prop error paracetamol-glucuronide in larva
 0 0.03749 0.06675  ; prop error paracetamol-sulphate in larva

$SIGMA  
 0 FIX  ; prop error paracetamol in blood
 0.09   ; prop error paracetamol-glucuronide in blood
 0.09   ; prop error paracetamol-sulphate in blood
 
 0.5    ; add error paracetamol in larva
 0.5    ; add error paracetamol-glucuronide in larva
 0 FIX  ; add error paracetamol-sulphate in larva
 
 0.05   ; add error paracetamol in blood
 0 FIX  ; add error paracetamol-glucuronide in blood
 0 FIX  ; add error paracetamol-sulphate in blood
 
 0 FIX  ; prop error paracetamol in medium
 0.09   ; prop error paracetamol-glucuronide in medium
 0.09   ; prop error paracetamol-sulphate in medium
 0.001  ; add error paracetamol in medium
 0 FIX  ; add error paracetamol-glucuronide in medium
 0.05   ; add error paracetamol-sulphate in medium
  
$ESTIMATION METHOD=1 MAXEVAL=8000 NOABORT PRINT=10 NSIG=3 SIGL=9 POSTHOC
$COVARIANCE PRINT=E
$TABLE ID TIME AMT DV EVID MDV CMT XEXP BQL XEXP2 AGE BLOOD IPRED
CWRES V2 V3 V4 K25 K12 K23 K24 K36 K47 T50 F5 F6 F7 ET1 ONEHEADER
NOPRINT FILE=tab_Zebrafish_paracetamol_metabolite
