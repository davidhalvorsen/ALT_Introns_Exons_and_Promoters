#!/usr/bin/env python
starting_length = 15000
division_loss = 50
current_length = starting_length
current_division_number = 1

while current_length > 0:
    if current_length > 5000:
        print "Telomere Length is " + str(current_length) + " bp" + " after " + str(current_division_number) + " divisions."
    elif current_length > 3000:
        print "p53-dependent arrest triggered @ " + str(current_length) + " bp" + " after " + str(current_division_number) + " divisions."
        print "BUT, this cell has transforming mutations! It continues to divide!"
    elif current_length > 0:
        print "Genomic instability is causing new mutations at " + str(current_length) + " bp" + " after " + str(current_division_number) + " divisions."
        print "Death is inevitable, UNLESS telomere length is stabilized!"
    else:
        print "David is a bad programmer, lol"
    current_length = current_length - 50
    current_division_number += 1
print "The telomere is gone!"
