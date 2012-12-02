#!/usr/bin/python

import sys
import os
import string
import shlex
import subprocess

prices = {}

def main():
  spot_prices = subprocess.check_output(['ec2-describe-spot-price-history', '-d Linux/UNIX (Amazon VPC) ','-d Linux/UNIX'])
  for line in spot_prices.split("\n"):    
    # sys.stdout.write("."+line + ".\n")
    if len(line) > 10: 
      (spot, price, timestamp, instance_type, product_description, zone) = line.split('\t')
      price = float(price)
      if instance_type not in prices.keys():
        prices[instance_type] = {}
        prices[instance_type]["MAX"] = price
        prices[instance_type]["MIN"] = price
        prices[instance_type]["TOTAL"] = price
        prices[instance_type]["COUNT"] = 1
      else:
        if prices[instance_type]["MAX"] < price: prices[instance_type]["MAX"] = price 
        if prices[instance_type]["MIN"] > price: prices[instance_type]["MIN"] = price
        prices[instance_type]["TOTAL"] += price
        prices[instance_type]["COUNT"] += 1

  sys.stdout.write("\t".join(["TYPE","MIN","MEAN","MAX\n"]))
  
  for instance_type in sorted(prices.keys()):
    sys.stdout.write("\t".join([instance_type,
                               str(prices[instance_type]["MIN"]),
                               str(prices[instance_type]["TOTAL"] / prices[instance_type]["COUNT"]),
                               str(prices[instance_type]["MAX"])+"\n"]))

if __name__ == "__main__":
  main()


