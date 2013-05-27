require "#{File.dirname(__FILE__)}/seq_lib.rb"
REF  = open("#{File.dirname(__FILE__)}/simulated.20130527.ref").read.gsub(/\n/,"")
DNA = ["A", "T", "G", "C"]
READ_SIZE=100
CHR = "chr1"
POS = -1
type = ARGV[0] #"nashi", "sm", "hap4"
depth = ARGV[1].to_i
tp = ARGV[2].to_f #tumor_purity
@out_f = ARGV[3]
BASE_QUAL = ARGV[4].to_i
@a = ARGV[5].to_i
@b = ARGV[6].to_i
@k = ARGV[7].to_i
@fn = [0.5, 0.5, 0.0, 0.0] # h1, h2, h3, h4
sum_f = tp * (@a + @b - 2) + 2
h1_f = ((1 - tp) + tp * (@a - @k)) / sum_f
h2_f = ((1 - tp) + tp * @b) / sum_f
h3_f = tp * @k / sum_f
warn ARGV

if type == "nashi"
  @ft = [0.5, 0.5, 0.0, 0.0]
elsif type == "sm"
  @ft = [h1_f, h2_f, h3_f, 0.0]
  #@ft = [(1.0-tp)/2.0, 1.0/2.0, tp/2.0, 0.0]
else
  @ft = [h1_f, h2_f, h3_f/2.0, h3_f/2.0]
  #@ft = [(1.0-tp)/2.0, 1.0/2.0, tp/4.0, tp/4.0]
end

p @fn
p @ft

def flip(pos)
  i = DNA.index(REF[pos])
  return i == DNA.length-1 ? DNA[0] : DNA[i+1]
end

snp_pos = 5000 + rand(5000)
@h1 = Haplotype.new
snp=Variant.new(snp_pos+1, "#{REF[snp_pos]}=>#{flip(snp_pos)}")
@h2 = Haplotype.new([snp])
sm_pos = snp_pos + rand(100) - 50
sm_pos += ((sm_pos - snp_pos).abs < 5) ? 20 : 0
sm=Variant.new(sm_pos+1,"#{REF[sm_pos]}=>#{flip(sm_pos)}")
@h3 = Haplotype.new([sm])
@h4 = Haplotype.new([sm, snp])
`echo "#{POS + sm_pos + 1}" > ans` if type == "sm"

if rand() > 0.5
  @haps = [@h1, @h2, @h3, @h4]
else
  @haps = [@h2, @h1, @h4, @h3]
end
p @haps

tumor = File.open(@out_f + "tumor.sam", "w")
normal = File.open(@out_f + "normal.sam", "w")
tumor.puts Haplotype::HEADER
normal.puts Haplotype::HEADER

if sm_pos < snp_pos
  l_pos = sm_pos;r_pos = snp_pos;
else
  l_pos = snp_pos;r_pos = sm_pos;
end

@normal_reads = []
index = 0
(depth*(200 + r_pos-l_pos).to_f/READ_SIZE).to_i.times do
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
  index+=1
  pair = @haps[i].gen_paired_reads(l_pos, r_pos)
  normal.puts Haplotype.pair_to_sam(index, pair)
end

index = 0
@tumor_reads = []
(depth*((200 + r_pos-l_pos).to_f/READ_SIZE)).to_i.times do # TODO: CNV should be considered
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
  index+=1
  pair = @haps[i].gen_paired_reads(l_pos, r_pos)
  tumor.puts Haplotype.pair_to_sam(index, pair)
end

