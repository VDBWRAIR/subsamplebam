def asbin(n):
        """converted a number to its binary rep (padded with 0's)"""
        return str(bin(n))[2:].zfill(17)

print "value\thex\tbinary"
for pow in range(17):
        val = 2 ** pow
        print "%-5d\t%-4x\t%s" % (val, val, asbin(val))

            # set all flags 
all_ones = reduce(lambda x, y: x | 2**y, range(17), 1)
print "\nall flags set:", asbin(all_ones)

