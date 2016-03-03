import math
import random
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import islice
import time
import pylab as P
# import plotly.plotly as py

"""
@Author: Tony Poerio
@email:  adp59@pitt.edu
University of Pittsburgh
Spring 2016
CS1538 - Simulation
Assignment #3 - Random Variate Generation

This python file will generate a series of random numbers from the standard normal distribution.
It will then compare them on the quality of the values they generated and the efficiency of the methods.
"""



####################
### CONTROL FLOW ###
####################
def main():
    print "UNIVERISTY OF PITTSBURGH - SPRING 2016:: CS1538, Assignment #3"
    print "--------------------------------------------------------------"

    # Get input, strip newlines.
    test_selection = ""

    while (test_selection != "q" ):
        select_test()
        test_selection = raw_input("Selection > ").strip()
        if test_selection == "q":
            exit()

        select_number_of_observations()
        number_observations = raw_input("Selection > ").strip()
        number_observations = int(number_observations)

        # If use selects python rand function,
        # create output file, and run the battery of tests
        if int(test_selection) == 1:
            inverse_transform( number_observations )
            time_experiment(test_selection, number_observations)
            chi_sq_result( int(test_selection), number_observations)

        # If use selects LCG function,
        # create output file, and run the battery of tests
        elif int(test_selection) == 2:
            accept_or_reject( number_observations )
            time_experiment(test_selection, number_observations)
            chi_sq_result( int(test_selection), number_observations)

        # If user selects LCG with RANDU settings
        # create output file, and run the battery of tests
        elif int(test_selection) == 3:
            polar_coords( number_observations )
            time_experiment(test_selection, number_observations)
            chi_sq_result( int(test_selection), number_observations)

        else:
            print "Please select a number from 1 to 3."

        # create_histogram( int(test_selection) )
        # TESTS
        # inverse_transform( 100 )
        # accept_or_reject( 100 )
        # polar_coords( 100 )


    return


##########################
### GENERATION METHODS ###
##########################


def inverse_transform( num_iterations ):
    """
    Yields numbers from about -4.0 to 4.0, with some larger outliers.
    These are Z-Scores.
    :param num_iterations: How many values to create (z-scores)
    :return: void. Writes values to a file named 'inverse_transform.txt'.
             File will be stored in the current directory.
    """
    # First, get CDF of the distro we need. Here that's Standard Normal.
    # Since Standard Normal CDF has no closed form, we'll use Bowling's 2009 approximation:
    #      Pr(Z <= z) = 1 / [1 + e^(-1.702z)
    # inverse obtained at:  http://www.wolframalpha.com/input/?i=solve+for+z:++R+%3D+1%2F(1%2Be%5E(1.072*z))
    e = math.e   # Get variable for "e"

    # counter for how many iterations we've run
    counter = 0

    # Open a file for output
    outFile = open("inverse_transform.txt", "wb")

    #Perfom number of iterations requested by user
    while counter < num_iterations:
        # Store value of each iteration
        R = random.random()   # Get a random number
        # cdf = 1 / (1 + e.__pow__(-1.702*z_score))  # get an approx for our cdf function
        x_value = (125.0/134.0)*math.log( (1.0-R)/R, e )    # get an inverse transformed function
                                                            # for z, in terms of r

        writeValue = str(x_value)
        # write to output file
        outFile.write(writeValue + "\n")
        print "num: " + " " + str(counter) +":: " + str(x_value)

        counter = counter+1

    outFile.close()
    print "Successfully stored " + str(num_iterations) + " random numbers in file named: 'inverse_transform.txt'."
    return



def accept_or_reject( num_iterations ):
    # RAND =  a draw from the U[0,1)
    RAND = random.random
    accept = False
    reject = False
    rejections = 0
    Z = 0      # The z-value we are generating


    # f(x) = the standard normal pdf function
    # f_of_x = (2.0 / math.sqrt(2.0 * math.pi)) * math.e.__pow__(- (RAND ** 2) / 2.0)

    # g(x) = the exponential density function with mean 1, that is: (lambda=1)
    # g_of_x = math.e.__pow__(-RAND)

    # c = max{f(x)/g(x)}
    # Formula from Sheldon Ross's Simulation: 5th addition
    c = math.sqrt( (2.0 * math.e) / math.pi )

    counter = 0

    # Open a file for output
    outFile = open("accept_or_reject.txt", "wb")

    #Perfom number of iterations requested by user
    while counter < num_iterations:

        # PROCEDURE, From ROSS: Simulation (5th Edition) Page 78
        # Step 1:  Generate Y1, an exponential random variable with rate 1
        Y1 = gen_exponential_distro_rand_variable()
        # Step 2:  Generate Y2, an exponential random variable with rate 2
        Y2 = gen_exponential_distro_rand_variable()
        # Step 3:  If Y2 - (Y1 - 1)^2/2 > 0, set Y = Y2 - (Y1 - 1)^2/2, and go to Step 4 (accept)
        #          Otherwise, go to Step 1 (reject)
        subtraction_value = ( math.pow( ( Y1 - 1 ), 2 ) ) / 2
        critical_value = Y2 - subtraction_value
        if critical_value > 0:
            accept = True
        else:
            reject = True

        # Step 4:  Generate a random number, U, and set:
        #          Z = Y1 if U <= 1/2
        #          Z = Y2 if U >- 1/2
        if accept == True:
            U = random.random()
            if (U > 0.5):
                Z = Y1
            if (U <= 0.5):
                Z = -1.0 * Y1

            writeValue = str(Z)
            # write to output file
            outFile.write(writeValue + "\n")
            print "num: " + " " + str(counter) +":: " + str(Z)
            counter += 1

        # Else, increment our rejection count
        else:
            rejections += 1


        # Reset boolean values
        accept = False
        reject = False


    outFile.close()
    print "Found this many rejections: " + str(rejections)
    print "Successfully stored " + str(num_iterations) + " random numbers in file named: 'accept_or_reject.txt'."

    return



def polar_coords( num_iterations ):
    reject = False
    rejections = 0
    writeValue1 = ""
    writeValue2 = ""
    X = 0.0
    Y = 0.0

    # counter for how many iterations we've run
    counter = 0

    # Open a file for output
    outFile = open("polar_coords.txt", "wb")

    #Perfom number of iterations requested by user
    while counter < num_iterations/2:
        # Process from Sheldon Ross' Simulations 5th edition, page 83
        # Step 1:  Generate random numbers U1 and U2
        U1 = random.random()
        U2 = random.random()
        # Step 2:  Set V1 = 2U1 - 1, V2 = 2U2 - 1, S = V1^2 + V2^2
        V1 = 2.0*U1 - 1.0
        V2 = 2.0*U2 - 1.0
        S = (V1 ** 2) + (V2 ** 2)

        if S > 1.0:
            reject = True
            rejections += 1
        if reject == False:
            X = math.sqrt( (-2.0 * math.log(S))/S )*V1
            Y = math.sqrt( (-2.0 * math.log(S))/S )*V2
            writeValue1 = str(X)
            writeValue2 = str(Y)
        # write to output file
        outFile.write(writeValue1 + "\n")
        outFile.write(writeValue2 + "\n")
        print "num: " + " " + str(counter) +":: " + str(X)
        print "num: " + " " + str(counter) +":: " + str(Y)

        counter = counter+1
        reject = False

    outFile.close()

    print "Found this many rejections: " + str(rejections)
    print "Successfully stored " + str(num_iterations) + " random numbers in file named: 'polar_coords.txt'."

    return

##############################
### STAT TESTS and OUTPUTS ###
##############################
def create_histogram( test_selection ):
    data_file = ""
    filename = ""
    title = ""
    if test_selection == 1:
        data_file = "inverse_transform.txt"
        filename = "inverse_transform"
        title = "INVERSE TRANSFORM"
    if test_selection == 2:
        data_file = "accept_or_reject.txt"
        filename = "accept_or_reject"
        title = "ACCEPT/REJECT"
    if test_selection == 3:
        data_file = "polar_coords.txt"
        filename = "polar_coords"
        title =  "POLAR COORDS"

    # Create histogram
    data_set = []  # dataset histogram


    # grabs first 100 files, as strings with newline endpoints
    with open( data_file, "r" ) as f:
        all_vals_as_STRINGS = list(islice(f, 10000))

    # transform all values to floats
    for val in all_vals_as_STRINGS:
        val = float(val)
        data_set.append(val)

    print ""
    #with open(data_file, "r") as f:
    #    # data points is a list containing all numbers we've read in.
    #    data_set = f.readlines()

    #for datapoint in data_set:
    #    datapoint = int(datapoint)

    # the histogram of the data with histtype='step'
    n, bins, patches = P.hist(data_set, 10, (-5, 5), normed=1, histtype='stepfilled')
    P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    P.title( title )
    P.xlabel( "Value Range" )
    P.ylabel( "Percentage" )
    #P.figure()
    P.show()
    P.savefig()

    #plt.hist(data_set)
    #plt.title(test_selection)
    #plt.xlabel("Value")
    #plt.ylabel("Num Observations")
    #fig = plt.gcf()
    # plot_url = py.plot_mpl(fig, filename=filename)

    print "created histogram for " + data_file + " stored as " + filename + "."
    print ""
    print ""



def chi_square_uniformity_test( data_set, confidence_level, num_samples ):
    """
    Null hypothesis:  Our numbers distributed uniformly on the interval [0, 1).

    This function uses the chi-square test for uniformity to determine whether our numbers
    are uniformly distributed on the interval [0,1).

    Formula is: "sum[ (observed-val - expected-val)^2 / expected val ], from 0 to num_samples"
    This gives us a number which we can test against a chi-square value table.

    Also need to know, degrees of freedom:  df=num_samples-1
    :param data_set: the data_set, must be a dictionary with 10 intervals.
                     Use return value from  @divide_RNG_data_into_10_equal_subdivisions_and_count
    :param confidence_level: confidence level we are testing at
    :param num_samples: number of data points
    :return: A chi-squared value
    """
    # This is our test statistic, this will be an accumulated value, as we loop through the data set
    chi_sq_value = 0.0
    degrees_of_freedom = num_samples - 1

    # We're doing 10 equal subdivisions, so need to divide our number samples by 10,
    # Assuming uniform distribution, to get an expected value. All values should be same
    # If our distro is actually uniform.
    expected_val = 0.0


    # Loop through a dictionary and get every count
    # The observed value is going to be our count at each key, and then we can do chi-square
    for observed_val in data_set:
        # print "Observed value is: " + observed_val
        difference = expected_val - observed_val
        chi_sq_value += ( math.pow((expected_val - observed_val), 2)/1 )

    # Coming out of this loop, we'll have a chi-squared test statistic
    # Now we just need to do a lookup to see if it's valid
    return chi_sq_value

def chi_sq_significance_test( chi_sq, signif_level):
    """
    Performs a significance test for df=10000, based on values calculated at:
    https://www.swogstat.org/stat/public/chisq_calculator.htm
    :param chi_sq:  Chi-sq value to test
    :param signif_level: Level of significance we are testing: 0.80, 0.90, or 0.95
    :return: message stating whether we accept or reject null
    """
    result = "FAIL TO REJECT null hypothesis"
    crit_value = 0.0
    if signif_level == 0.8:
        crit_value = 10118.8246
    elif signif_level == 0.90:
        crit_value = 10181.6616
    elif signif_level == 0.95:
        crit_value = 10233.7489
    else:
        print "**Invalid Significance Level for Chi Sq***"

    if chi_sq > crit_value:
        result = "REJECT null hypothesis"

    print "Print Significance Level: " + str(signif_level)
    print "Chi Sq: " + str(chi_sq)
    print "Crit Value: " + str(crit_value)
    print "Result is: " + result
    print "...................................."

    return result

def divide_RNG_data_into_10_equal_subdivisions_and_count( data_file ):
    """
    Takes a path to a data file in the current directory.
    Returns a dictionary with keys 1-10, values=num instances in each of
    10 equal intervals from range: [0, 1).
    The function counts how many data points are in each interval, and gives us
    a dictionary so we can manipulate this data more easily, based on count by index.

    :param data_file: Must be in current directory. Pass in the string name.
    :return: A dictionary with counts of how many occurrences our data had for each
    of 10 equal intervals between [0, 1). (Divided into 10ths)
    """
    # For each of our uniformity tests, need to divide our data points in 10 equal subdivisions
    subdivisions = {  "1":  0,
                      "2":  0,
                      "3":  0,
                      "4":  0,
                      "5":  0,
                      "6":  0,
                      "7":  0,
                      "8":  0,
                      "9":  0,
                      "10": 0   }
    with open(data_file, "r") as f:
        # data points is a list containing all numbers we've read in.
        data_points = f.readlines()

    # Loop through our data points and count number of data points in each subdivision
    # Divide by tenths, from 0.0 to 1.0.
    for num in data_points:
        num = float(num)
        if num < -4.0:
            subdivisions["1"] += 1
        elif num < -3.0:
            subdivisions["2"] += 1
        elif num < -2.0:
            subdivisions["3"] += 1
        elif num < -1.0:
            subdivisions["4"] += 1
        elif num < 0.0:
            subdivisions["5"] += 1
        elif num < 1.0:
            subdivisions["6"] += 1
        elif num < 2.0:
            subdivisions["7"] += 1
        elif num < 3.0:
            subdivisions["8"] += 1
        elif num < 4.0:
            subdivisions["9"] += 1
        elif num < 5.0:
            subdivisions["10"] += 1

    return subdivisions


######################
### HELPER METHODS ###
######################

def time_experiment( generation_method, num_iterations ):
    begin_time = time.time()
    if generation_method == 1:
        inverse_transform( num_iterations )
    elif generation_method == 2:
        accept_or_reject( num_iterations )
    elif generation_method == 3:
        polar_coords( num_iterations )
    elif generation_method < 1 or generation_method > 3:
        print "Invalid method entered. Please select a number from 1 to 3."

    end_time = time.time()

    total_time = begin_time - end_time

    print "Total time spent: " + str(total_time)
    print ""
    print ""
    return total_time


def gen_normal_pdf_rand_variable():
    return (2 / math.sqrt(2 * math.pi)) * math.e.__pow__(- (random.random() ** 2) / 2)

def gen_exponential_distro_rand_variable():
    # return math.e.__pow__( -1 * random.random() )
    return -( math.log( (  random.random() ) ) )

def chi_sq_result( test_selection, num_samples ):
    data_file = ""
    filename = ""
    title = ""
    if test_selection == 1:
        data_file = "inverse_transform.txt"
        filename = "inverse_transform"
        title = "INVERSE TRANSFORM"
    if test_selection == 2:
        data_file = "accept_or_reject.txt"
        filename = "accept_or_reject"
        title = "ACCEPT/REJECT"
    if test_selection == 3:
        data_file = "polar_coords.txt"
        filename = "polar_coords"
        title =  "POLAR COORDS"

    # Create histogram
    data_set = []  # dataset histogram

    # grabs first 100 files, as strings with newline endpoints
    with open( data_file, "r" ) as f:
        all_vals_as_STRINGS = list(islice(f, 10000))

    # transform all values to floats
    for val in all_vals_as_STRINGS:
        val = float(val)
        data_set.append(val)

    chi_sq = chi_square_uniformity_test( data_set, 0.05, num_samples )
    print "CHI SQ RESULT IS: " + str(chi_sq)
    chi_sq_significance_test(chi_sq, 0.95)


def select_test():
    """
    Command line prompt for selecting a test
    :return: void - prints a prompt to command line
    """
    print "Please select a method for generating random variates for the Standard Normal Distribution: "
    print " 1. Inverse Transform Method "
    print " 2. Accept/Reject Method "
    print " 3. Polar-Coordinate Method "
    print ""
    print "      (or type 'q' to quit)"
    print ""

def select_number_of_observations():
    """
    Command line prompt to select the number of observations for a given test
    :return: void - prints a prompt to command line
    """
    print "How many observations should we perform?"


def run_test_suite( test_selection, number_observations ):
    """
    Runs all of our test suites and prints output to the screen
    :param test_selection:  an int - 1,2, or 3. Corresponds to test selected.
    :param number_observations: the number of data points to test
    :return: void - prints to command line
    """
    begin_time = time.time()
    input_file = ""
    test_name = ""
    test_selection = int(test_selection)
    if test_selection == 1:
        input_file = "py_random_output.txt"
        test_name = "PYTHON BUILT-IN RAND"

    elif test_selection == 2:
        input_file = "lgc_output.txt"
        test_name = "LINEAR CONGRUENTIAL GENERATOR"

    elif test_selection == 3:
        input_file = "lgc_RANDU_output.txt"
        test_name = "LGC with RANDU initial settings"
    else:
        print "Invalid input. Please try again."

    print ""
    print ""
    print "TEST SUITE FOR:  %s " % (test_name)
    print "======================================"


###################
### ENTRY POINT ###
###################

if __name__ == "__main__":
    main()

