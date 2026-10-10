/ kdb+tick schema: time and sym must be the first two columns
quote:([]time:`timespan$();sym:`g#`symbol$();bid:`float$();ask:`float$();bsize:`long$();asize:`long$())
