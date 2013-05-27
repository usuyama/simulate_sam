require "#{File.dirname(__FILE__)}/seq_lib.rb"
REF  = open("#{File.dirname(__FILE__)}/simulated.20130526.ref").read.gsub(/\n/,"")
DNA = ["A", "T", "G", "C"]
READ_SIZE=100
type = ARGV[0] #"nashi", "hetero_snp", "sm"
depth = ARGV[1].to_i
tp = ARGV[2].to_f #tumor_purity
#del_len = ARGV[2].to_i
@out_f = ARGV[3] || "./"
@fn = [0.5, 0.5, 0.0, 0.0] # h1, h2, h3, h4
@ft = [(1.0-tp)/2.0, 1.0/2.0, tp/2.0, 0.0]

p @fn
p @ft

def flip(pos)
  i = DNA.index(REF[pos])
  return i == DNA.length-1 ? DNA[0] : DNA[i+1]
end
snp_pos = 1400 + rand(200)
@h1 = Haplotype.new
p Variant.new(snp_pos+1, "#{REF[snp_pos]}=>#{flip(snp_pos)}")
@h2 = Haplotype.new([Variant.new(snp_pos+1, "#{REF[snp_pos]}=>#{flip(snp_pos)}")])
sm_pos = snp_pos + rand(120) - 60
sm_pos += (sm_pos == snp_pos) ? 10 : 0
p Variant.new(sm_pos+1,"#{REF[sm_pos]}=>#{flip(sm_pos)}")
@h3 = Haplotype.new([Variant.new(sm_pos+1,"#{REF[sm_pos]}=>#{flip(sm_pos)}")])
@h4 = Haplotype.new([Variant.new(snp_pos+1, "#{REF[snp_pos]}=>#{flip(snp_pos)}"), Variant.new(sm_pos+1,"#{REF[sm_pos]}=>#{flip(sm_pos)}")])
`echo "#{10191429 + 40 + sm_pos + 1}" > ans` if type == "sm"

if type == "nashi"
  @haps = [ @h1, @h2]
elsif type == "hetero_snp"
  @haps = [ @h2, @h3 ]
else
  if rand() > 0.5
    @haps = [	@h1, @h2, @h3, @h4]
  else
    @haps = [	@h2, @h1, @h4, @h3]
  end
end
p @haps

tumor1 = File.open(@out_f + "tumor1", "w")
tumor2 = File.open(@out_f + "tumor2", "w")
normal1 = File.open(@out_f + "normal1", "w")
normal2 = File.open(@out_f + "normal2", "w")

@normal_reads = []
(depth*3).times do
  r = rand()
  i = 0
  th = 0.0
  0.upto(@haps.size()-1) do
    th += @fn[i]
    if th >= r
      break
    else
      i += 1
    end
  end
  pair = @haps[i].gen_paired_fastq(snp_pos)
  normal1.puts(pair[0])
  normal2.puts(pair[1])
end

@tumor_reads = []
(depth*3).times do
  r = rand()
  i = 0
  th = 0.0
  0.upto(@haps.size()-1) do
    th += @ft[i]
    if th >= r
      break
    else
      i += 1
    end
  end
  pair = @haps[i].gen_paired_fastq(snp_pos)
  tumor1.puts(pair[0])
  tumor2.puts(pair[1])
end

