#!/usr/bin/env python
starting_length = 15000
division_loss = 50
current_length = starting_length
current_division_number = 1

while current_length > 0:
    print "Telomere Length is " + str(current_length) + " bp" + " after " + str(current_division_number) + " divisions."
    current_length = current_length - 50
    current_division_number += 1
print "The telomere is gone!"
