; ------------------------------------------------------------------------------------------
;
; An example of Genetic Expression Programming in NetLogo
;
; Author:   jihe.Gao@jiejiaotech.com
;
; ------------------------------------------------------------------------------------------
extensions [csv table time]

breed [ genes gene ]
breed [ codeturtles codeturtle ]

globals [
  gp_weighted_gene_list
  gp_syntaxlist
  gp_generation
  avgfitness
  bestfitness
  worstfitness
  gp_best_fitness_this_gen
  the_watched_one
]


genes-own
[
  c_type
  input_var_list
  output_type  ; void (command), reporter, boolean
  odds
]

codeturtles-own
[
  my_code
  my_tree
  gp_raw_fitness
  birthtime
  result
]





to setup
  clear-all
  reset-ticks
  gp-setup
end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                   START OF GP LIBRARY CODE                         ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;      gpuser-setup-syntax    (determines the ingredients for randomlly generated code trees)
;;      gpuser-run-codeturtles  (determines what happens to the codeturtles during a generation)
;;      gpuser-initialize-codeturtle  (run when each codeturtle that is created, sets initial variables, etc)
;;      gpuser-calculate-raw-fitness   (reports some measure analogous to "distance to goal", where 0 is perfection)
;;      gpuser-calculate-adjusted-fitness (reports a value between 0 and 1, where 1 is perfection)
;;
;; In addition, your setup should call "gp-setup" and your go should call "gp-episode"
;;
;; Notes
;;
;;  Structure of a tree node:
;;    ["NODE_VALUE"  RTYPE_*  TYPE_*  CHILD1  CHILD2  ... ]
;;


to-report gp-target-answer
  report ticks * 3
end



to gp-setup
  gp-setup-genes
  gp-setup-codeturtles
end


to gp-setup-genes
  ;set gp_syntaxlist csv:from-file "genes.csv"
  set gp_syntaxlist init-gp-syntaxlist
  foreach gp_syntaxlist [ [?]->
    if not any? turtles with [label = first ?][
      create-genes 1 [
        ;["gp_if" "TYPE_COMMAND" "RTYPE_BOOLEAN" "RTYPE_COMMAND" 10]
        set label first ?
        set input_var_list sublist ? 2 (length ? - 1)
        set output_type item 1 ?
        set odds last ?
        move-to one-of patches with [pycor = min-pycor]
      ]
    ]
  ]
  set gp_weighted_gene_list (reduce [[?1 ?2]-> (sentence ?1 ?2)] (map [[?]-> n-values ([odds] of ?) [[label] of ?] ] sort genes))
end


to gp-setup-codeturtles
  create-codeturtles population-size [ gp-initialize-codeturtle [0]]
end


;; whenever codeturtles are created (at the beginning of each generation) this procedure is run
;; gen_config: (list [my_code] of winner [my_code] of partner )
to gp-initialize-codeturtle [gene_lists]
  ifelse gene_lists = [0]
  [
    set my_code (gp-weighted-random-genes initial-chromosome-length)
  ]
  [
    set my_code (gp-cross-chromosomes gene_lists)
  ]
  carefully [
    set my_tree gp-compile-tree (sublist my_code (1 + position "Q" my_code) length my_code)
  ][set my_tree "0"]
  set birthtime ticks
  set size 2
  set shape "circle"
  set gp_raw_fitness 100  ; large number (very unfit)
  set xcor random-xcor ; find the start patch
end


;; reports a random entry from a syntaxlist, where the chance of an entry being chosen
;;  is weighted by the number that is the last entry of the element
to-report gp-weighted-random-genes [ length_gene ]
  report fput "Q" n-of length_gene gp_weighted_gene_list
end


to-report gp-cross-chromosomes [DNA_lists]
  let DNA0 first DNA_lists
  let DNA1 last  DNA_lists
  ifelse length DNA0 < length DNA1
  [ set DNA1 sublist DNA1 0 (1 + length DNA0) ]
  [ if length DNA0 > length DNA1 [
    set DNA1 (sentence DNA1 n-values (length DNA0 - length DNA1) ["0"]) ]
  ]
  let dna []
  (foreach DNA0 DNA1 [ [d1 d2]->
    ifelse (random 100 < crossover-chance)
    [ set dna lput d2 dna ]
    [
      ifelse (random 100 < mutate-chance)
      [ set dna lput (one-of gp_weighted_gene_list) dna ]
      [ set dna lput d1 dna ]]
  ])
  report dna
end


;;
to-report gp-compile-tree [this_code]
  ; this_code : [ "gp_+" "3" "8" "32" "0" ]

  let thegene one-of turtles with [label = first this_code]
  set this_code but-first this_code
  let tree reduce [[?1 ?2]-> (word ?1 " " ?2)] [(sentence label input_var_list)] of thegene

  while [not empty? this_code and member? "R" tree][
    set thegene one-of turtles with [label = first this_code]
    set this_code but-first this_code
    let index position "R" tree
    set tree replace-item index tree (insert-node-string thegene)
  ]
  report tree
end

to-report insert-node-string [thegene]
  report reduce [[?1 ?2]-> (word " " ?1 " " ?2 " ")] [(sentence label map [[?]-> (word "(" ? ")")] input_var_list)] of thegene
end







;; This is the main "go" procedure for the gp-library.
to gp-episode
  gp-create-next-generation ; breed_option = "wta" "equal"
  gp-run-codeturtles
  gp-update-fitness
  if move? [ask codeturtles [set ycor min-pycor + gp-fitness * world-height - 1]]
  gp-plot
  tick
end



;; codeturtles procedure
;; Think of this reporter as "how far am I from my goal?"
;; Thus, a perfect codeturtle has a raw fitness of 0
;;       a bad codeturtle has an arbitrarily large raw fitness

to gp-run-codeturtles
  ask codeturtles [carefully [set result (runresult my_tree)][die]]
end



to gp-update-fitness
  ask codeturtles [
    set gp_raw_fitness gp-calculate-fitness
  ]
end

; user-defined fitness function ( cost function )
to-report gp-calculate-fitness
  report abs (gp-target-answer - result)
end

to-report gp-fitness
  report 1 / (1 + gp_raw_fitness)
end



to gp-create-next-generation

  if breed_option = "wta" ; winner takes all
  [
    let winner min-one-of codeturtles [gp_raw_fitness]
    repeat population-size - count turtles [
      let partner one-of codeturtles
      gp-breed-new-codeturtle (list gp_generation [my_code] of winner [my_code] of partner ) ]
  ]
  if breed_option = "equal"
  [
    repeat population-size - count turtles [
      let couple n-of 2 codeturtles with [gp-fitness >= avgfitness]
      if couple = nobody [set couple n-of 2 codeturtles]
      gp-breed-new-codeturtle (list [my_code] of first sort couple [my_code] of last sort couple) ]
  ]

if count codeturtles > n-die-out-each-tick [ ask n-of n-die-out-each-tick codeturtles [ die ] ]
end


; gen_config: (list [my_code] of turtle1 [my_code] of turtle2 )
to gp-breed-new-codeturtle [gene_lists]
  create-codeturtles 1 [ gp-initialize-codeturtle gene_lists ]
end







;; plot the fitnesses of the current generation, and update the gp-best-fitness-this-gen variable
to gp-plot
  if any? codeturtles [
    set bestfitness max [ gp-fitness ] of codeturtles
    set avgfitness mean [ gp-fitness ] of codeturtles
    set worstfitness min [ gp-fitness ] of codeturtles
    set gp_best_fitness_this_gen bestfitness

    set-current-plot "Fitness Plot"
    set-current-plot-pen "avg"
    plot avgfitness
    set-current-plot-pen "best"
    plot bestfitness
    set-current-plot-pen "worst"
    plot worstfitness
  ]
end



;; displays the code of the best-this-gen codeturtle in the output box
to gp-showbest
  clear-output
  ask one-of (codeturtles with-min [ gp_raw_fitness ])
  [ print (word self " tree: " my_tree " result: " result )]
end

to gp-show-all
  ;clear-output
  foreach sort-by [[?1 ?2]-> [gp_raw_fitness] of ?1 < [gp_raw_fitness] of ?2] codeturtles
  [
    [?]-> ask ? [
      show (word self " code: " my_code " tree: " my_tree " result: " result )
    ]
  ]
end


to select-one
  if mouse-down? [
    set the_watched_one closest-codeturtle mouse-xcor mouse-ycor turtles with [breed != links]
    ask links [die]
    ask turtles with [shape = "square"][die]
    ask codeturtles [set label ""]
    watch the_watched_one
    ask the_watched_one [
      set label my_tree
      setxy mouse-xcor mouse-ycor
    ]
    display
  ]
end

to-report closest-codeturtle [x y agent-set]  ; Return closest agent to x, y
 report min-one-of codeturtles [distancexy x y]
end
to reset-codeturtles

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                    END OF GP LIBRARY CODE                          ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; These are atom genes of procedures and reporters, that help form the DNA of the codeturtles.

to-report init-gp-syntaxlist
report [
  ["gp_+"      "TYPE_OPERATOR" "R" "R" 10]
  ["gp_-"      "TYPE_OPERATOR" "R" "R" 10]
  ["gp_*"      "TYPE_OPERATOR" "R" "R" 5]
  ["gp_safediv" "TYPE_REPORTER" "R" "R" 5]
  ["gp_random"  "TYPE_REPORTER" "R" 10]
  ["gp_ticks"   "TYPE_REPORTER" 5]
  ["0" "TYPE_REPORTER" 20]
  ["1" "TYPE_REPORTER" 20]
  ["2" "TYPE_REPORTER" 20]
  ["3" "TYPE_REPORTER" 20]
  ["4" "TYPE_REPORTER" 20]
  ["5" "TYPE_REPORTER" 20]
  ["8" "TYPE_REPORTER" 20]
  ["16" "TYPE_REPORTER" 20]
  ["32" "TYPE_REPORTER" 20]
]
end

to-report gp_ticks
  report ticks
end

to-report gp_ifelse-value [bool tree1 tree2]
  report ifelse-value bool [ runresult tree1 ][ runresult tree2 ]
end


to-report gp_+ [input1 input2]
  report input1 + input2
end

to-report gp_- [input1 input2]
  report input1 - input2
end

to-report gp_* [input1 input2]
  report input1 * input2
end

to-report gp_safediv [input1 input2]
  report ifelse-value (input2 != 0) [input1 / input2][999999]
end

to-report gp_random [input]
  report random input
end


to-report gp_>= [input1 input2]
  report input1 >= input2
end

to-report gp_<= [input1 input2]
  report input1 <= input2
end

to-report gp_> [input1 input2]
  report input1 >= input2
end

to-report gp_< [input1 input2]
  report input1 >= input2
end
@#$#@#$#@
GRAPHICS-WINDOW
10
15
565
571
-1
-1
8.97
1
15
1
1
1
0
1
1
1
-30
30
-30
30
0
0
1
ticks
30.0

BUTTON
582
22
684
55
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
582
184
1068
360
Fitness Plot
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"avg" 1.0 0 -2674135 true "" ""
"best" 1.0 0 -13345367 true "" ""
"worst" 1.0 0 -7500403 true "" ""

INPUTBOX
583
371
664
431
randomseed
2.0
1
0
Number

INPUTBOX
666
371
757
431
population-size
200.0
1
0
Number

INPUTBOX
584
436
789
496
initial-chromosome-length
15.0
1
0
Number

INPUTBOX
757
371
845
431
mutate-chance
10.0
1
0
Number

INPUTBOX
845
371
943
431
crossover-chance
30.0
1
0
Number

BUTTON
688
23
791
56
NIL
gp-episode
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
581
143
689
176
NIL
 gp-show-all
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
945
380
1068
425
breed_option
breed_option
"wta" "equal"
1

BUTTON
690
145
803
178
NIL
gp-showbest
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
792
436
944
496
n-die-out-each-tick
10.0
1
0
Number

SWITCH
581
104
684
137
move?
move?
1
1
-1000

MONITOR
849
93
942
138
NIL
count turtles
17
1
11

BUTTON
804
145
901
178
NIL
select-one
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
905
141
1025
186
NIL
the_watched_one
17
1
11

@#$#@#$#@
# Gene expression programming in NetLogo

## WHAT IS IT?

[Gene expression programming (GEP) ]([https://en.wikipedia.org/wiki/Gene_expression_programming]) is an evolutionary algorithm that creates computer programs or models. These computer programs are complex tree structures that learn and adapt by changing their sizes, shapes, and composition, much like a living organism. And like living organisms, the computer programs of GEP are also encoded in simple linear chromosomes of fixed length. Thus, GEP is a genotype-phenotype system, benefiting from a simple genome to keep and transmit the genetic information and a complex phenotype to explore the environment and adapt to it.

## HOW IT WORKS

You can find details in http://www.gene-expression-programming.com/webpapers/GEP.pdf


## HOW TO USE IT

0, (optional) modify "gp-target-answer" in code tab.
gp-target-answer is the target function which is the goal for the GEP to find out. 
By default, the target function is:

>to-report gp-target-answer
report ticks * 3
end

You can try 
> report 100

or
> report ticks * 2 + 100

etc.


1, click "setup" button (or type "setup" in command-center then hit ENTER)

2, click "gp-episode"

3, Observe the GEP process by watching "Fitness Plot", or hit "gp-show-all", "gp-showbest" or "select-one" to see code(s) and code tree(s). If the plot line "best" reached and stablized at 1, then there is some codeturtles successfully invented the target function. If the "best", "avg" is far below 1 at long run, you may change the option in  "breed_option" chooser to see if there is some improvement. 

(check model-vew.png for instance)

4, to add more code genes, first locate the comments in the code tab:
>;; These are atom genes of procedures and reporters, that help form the DNA of the codeturtles.

Then add your reporter function. 

Last, add the function name, type, "R"(parameter place holder) in the init-gp-syntaxlist reporter:

>to-report init-gp-syntaxlist
report [
>  ["gp_+"      "TYPE_OPERATOR" "R" "R" 10]
...

## CREDITS AND REFERENCES

https://github.com/jihegao/NetLogo_GEP

contact: jihe.gao(at)jiejiaotech.com
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
