# Gene expression programming in NetLogo

## WHAT IS IT?

[Gene expression programming (GEP) ]([https://en.wikipedia.org/wiki/Gene_expression_programming]) is an evolutionary algorithm that creates computer programs or models. These computer programs are complex tree structures that learn and adapt by changing their sizes, shapes, and composition, much like a living organism. And like living organisms, the computer programs of GEP are also encoded in simple linear chromosomes of fixed length. Thus, GEP is a genotype-phenotype system, benefiting from a simple genome to keep and transmit the genetic information and a complex phenotype to explore the environment and adapt to it.

## HOW IT WORKS

You can find details in http://www.gene-expression-programming.com/webpapers/GEP.pdf


## HOW TO USE IT

0, (optional) modify "gp-target-answer" in code tab, which is by default

>to-report gp-target-answer
report ticks * 3
end

gp-target-answer is the target function which is the goal for the GEP to find out.

1, click "setup" button (or run "setup" in command-center)
2, click "gp-episode"
3, Observe the GEP process by watching "Fitness Plot", or hit "gp-show-all", "gp-showbest" or "select-one" to see code(s) and code tree(s). If the plot line "best" reached and stablized at 1, then there is some codeturtles successfully invented the target function. If the "best", "avg" is far below 1 at long run, you may change the option in  "breed_option" chooser to see if there is some improvement. If still


## CREDITS AND REFERENCES

https://github.com/jihegao/NetLogo_GEP

contact: jihe.gao(at)jiejiaotech.com
