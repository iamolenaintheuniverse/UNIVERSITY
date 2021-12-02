from Pyro4 import expose
import datetime
import random


def modular_pow(base, exponent, mod):
    x = 1
    y = base
    while (exponent > 0):
        if (exponent % 2 == 1):
            x = (x * y) % mod
        y = (y * y) % mod
        exponent = exponent // 2

    return x % mod


def Jacobian(a, n):

    if (a == 0):
        return 0;# (0/n) = 0

    ans = 1
    if (a < 0):
        
        # (a/n) = (-a/n)*(-1/n)
        a = -a
        if (n % 4 == 3):
        
            # (-1/n) = -1 if n = 3 (mod 4)
            ans = -ans

    if (a == 1):
        return ans; # (1/n) = 1

    while (a):
        if (a < 0):
            
            # (a/n) = (-a/n)*(-1/n)
            a = -a
            if (n % 4 == 3):
                
                # (-1/n) = -1 if n = 3 (mod 4)
                ans = -ans

        while (a % 2 == 0):
            a = a // 2
            if (n % 8 == 3 or n % 8 == 5):
                ans = -ans

        a, n = n, a

        if (a % 4 == 3 and n % 4 == 3):
            ans = -ans
        a = a % n

        if (a > n // 2):
            a = a - n

    if (n == 1):
        return ans

    return 0

class Solver:

    def __init__(self, workers=None, input_file_name=None, output_file_name=None):
        self.input_file_name = input_file_name
        self.output_file_name = output_file_name
        self.workers = workers
        print("Inited")
    
    def solve(self):
        print("Job Started")
        print ('Number of Workers {0}'.format(len(self.workers)))

        (n, k) = self.read_input()
        NumIterPerWorker = k // len(self.workers)
        remainder = k % len(self.workers)
        lastWorkerNum = len(self.workers) - 1

        mapped = []

        for i in xrange(0, lastWorkerNum):
            print("map %d" % i)
            mapped.append(self.workers[i].mymap(str(n), NumIterPerWorker))
        mapped.append(self.workers[lastWorkerNum].solovoyStrassen(str(n), NumIterPerWorker + remainder))

        reduced = self.myreduce(mapped)

        self.write_output(reduced)

        print("Job Finished")
        
    @staticmethod
    @expose
    def mymap(n, k):
        print(n)
        print(k)
        a = int(n)
        prime = Solver.solovoyStrassen(n, k)

        return prime

    @staticmethod
    @expose
    def myreduce(mapped):
        print("reduce")
        output = []

        for isPrime in mapped:
            print("reduce loop")
            output.append(str(isPrime.value))
        print("reduce done")
        return output

    def read_input(self):
        f = open(self.input_file_name, 'r')
        n = int(f.readline())
        k = int(f.readline())
        f.close()

        return n, k

    def write_output(self, output):
        f = open(self.output_file_name, 'w')
        f.write('Solovay-Strassen Primality Test\n')
        if False in output:
             f.write('Result: n is COMPOSITE.\n')
        else:
            f.write('Result: n is Probably PRIME.\n')

        f.close()
        print("output done")

    @staticmethod
    @expose
    def solovoyStrassen(n, iterations):

        n = int(n)

        if (n < 2):
            return False
        if (n != 2 and n % 2 == 0):
            return False


        for i in xrange(iterations):
            a = random.randrange(n - 1) + 1
            jacobian = (n + Jacobian(a, n)) % n
            mod = modular_pow(a, (n - 1) / 2, n)

            if (jacobian == 0 or mod != jacobian):
                return False
        return True

if __name__ == '__main__':
   master = Solver([Solver(), Solver()], "CloudComputing/Solovay_Strassen/input2.txt", "CloudComputing/Solovay_Strassen/output.txt")
   master.solve()