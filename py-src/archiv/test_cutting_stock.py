#------------------------------------------------------------------------------
#> Test the cutting stock problem, applied to the batch system of HAWK.
#
#> @author OpenAI - ChatGPT
#> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
#> Date:    13.03.2023
#> LastMod: 13.03.2023
#
#> @brief:
#> Written by ChatGPT as an answer to the question
#> "Please write a python code to solve the cutting stock problem."
#>
#> Answer: Sure, here's an example Python code to solve the Cutting Stock 
#> Problem using the First Fit Decreasing (FFD) heuristic:
#> 
#> In this code, item_lengths is a list of the lengths of the items to be 
#> cut from a stock of length stock_length. The function cutting_stock 
#> returns the minimum number of stock pieces required to cut all the items 
#> using the First Fit Decreasing (FFD) heuristic.
#------------------------------------------------------------------------------
import math

def cutting_stock(item_lengths, stock_length):
    # sort the items in decreasing order
    item_lengths.sort(reverse=True)
    
    # initialize the number of bins used and the current bin capacity
    num_bins = 0
    current_bin_capacity = 0
    
    # iterate through the items
    for item_length in item_lengths:
        # try to fit the item in the current bin
        if current_bin_capacity >= item_length:
            current_bin_capacity -= item_length
        else:
            # if the item doesn't fit in the current bin, start a new bin
            num_bins += 1
            current_bin_capacity = stock_length - item_length
    
    return num_bins

# Example usage
item_lengths = [4, 3, 2, 1]
stock_length = 10
num_bins = cutting_stock(item_lengths, stock_length)
print("Number of bins used: ", num_bins)

