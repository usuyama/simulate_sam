#!/usr/bin/env ruby
# gen simulated fastq
# 51 chrs in 1 line
bases = ['A', 'T', 'G', 'C']
300.times do
  51.times do
    print bases[rand(4)] 
  end
  puts ''
end
