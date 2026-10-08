/ schema kdb+tick : time et sym obligatoirement en 2 premieres colonnes
quote:([]time:`timespan$();sym:`g#`symbol$();bid:`float$();ask:`float$();bsize:`long$();asize:`long$())
