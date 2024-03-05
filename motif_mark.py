#!/usr/bin/env python

import cairo
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="A proram to find motifs in a sequence")
    parser.add_argument("-f", "--fasta_file", help="to specify fasta input file", required=True)
    parser.add_argument("-m", "--motif_file", help="to specify the motif input file", required=True)
    parser.add_argument("-o", "--out_file", help="to specify the output image file")#required=True)
    return parser.parse_args()

args=get_args() 
fasta_file:str =args.fasta_file
motif_file:str =args.motif_file
png_file:str =args.out_file


def oneline_fasta(fasta_file):
    '''This function takes in a fast file and returns a dictionary with headers as keys and sequences as values'''
    fasta_dict = {}
    header = ''
    sequence = ''
    with open(fasta_file, 'r') as fa_file:
        for line in fa_file:
            line = line.strip()
            if line.startswith('>'):
                if header != '':
                    fasta_dict[header] = sequence
                    sequence = ''
                header = line
            else:
                sequence += line
        # Add the last entry to the dictionary
        if header != '':
            fasta_dict[header] = sequence
    return fasta_dict
                    
#create dictionary to find ambiguity
iupac = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'U':'[TU]',
         'W':'[ATU]','S':'[CG]', 'M':'[AC]', 'K':'[GTU]',
         'R':'[AG]', 'Y':'[CTU]','B':'[CGTU]','D':'[AGTU]',
          'H':'[ACTU]','V':'[ACG]', 'N':'[ACGTU]' }
#print(iupac)

class Gene:
    def __init__(self):
        '''creating the gene'''
        self.gene_name = None

    def get_name(self, key):
        #seq_dict = oneline_fasta(fasta_file)
        #for key in seq_dict:
        key_vals = key.split()
        self.gene_name = key_vals[0][1:]
        return self.gene_name

class Intron:
    def __init__(self):
        self.color = 0, 0, 0, 1 
        self.intron_length = None
    
    def get_intron_length(self,sequence):
        #fasta_seq = oneline_fasta(fasta_file)
        #for key, value in fasta_seq.items():
        #length = len()
        self.intron_length = len(sequence)
        return self.intron_length
    
class Exon:
    def __init__(self):
        self.color =  0, 0, 0
        self.exon_pos = None
    
    def get_exon_pos(self, sequence):
        #fasta_seq = oneline_fasta(fasta_file)
        self.exon_pos = []
        #for key, value in fasta_seq.items():
        for base in range(len(sequence)):
            if sequence[base].isupper():
                self.exon_pos.append(base)
        return self.exon_pos
            
        #return self.intron_length

class MotifClass:
    def __init__(self):
        '''motif definitaions'''
        self.color_pallette = [(144/255,238/255,144/255),(255/255,192/255,203/255),(128/255,0,128/255),(155/255,165/255,0),
                        (0,1,1),(128/255,128/255,128/255),(1,1,0),(128/255,128/255,0),(1,0,1),(1,20/255,147/255)]
        self.motifs = []
        self.match_positions = []
        self.regex_str = ''
        self.regex_list = []
        self.motif_colors = None
        self.actual_motif = None
        self.matchings = {}
        self.replace_regex = []
        self.motif_dictionary = {}
    
    def get_motifs(self, motif_file):
        with open(motif_file, 'r') as motifs_file:
            for line in motifs_file:
                line = line.strip().upper()
                self.motifs.append(line)
        return self.motifs
    
    def get_regex(self):
        #print(f'before isupper: {sequence}')
        #regex_str =''
        #sequence = motif_sequence.upper()
        # for motif in self.motifs:
        #     motif = motif.upper()
        #     print(motif)
        #     for base in motif:
        #         if base in iupac:
        #             #print('base in iupac: {base}')
        #             self.regex_str+=iupac[base]
        # return self.regex_str
        for motif in self.motifs:
            motif = motif.upper()
            replaced_motif = []

            for base in motif:
                if base in iupac:
                    replaced_motif.append(iupac[base])
                else:
                    replaced_motif.append(base)

            self.regex_str = ''.join(replaced_motif)
            self.regex_list.append(self.regex_str)
        return self.regex_list
    
   
    
    def get_replace_regex(self,sequence):
        #print(f'before isupper: {sequence}')
        #motif_sequence = self.motifs.upper()
        sequence = sequence.upper()
        for base in sequence:
            if base in iupac:
                self.replace_regex.append(iupac[base])
            else:
                self.replace_regex.append(base)

        return ''.join(self.replace_regex)
    
    def find_motifs(self, regex_dict, sequence):
        sequence = sequence.upper()
        pattern_parts = []
        new_dict ={}
        for key, value in regex_dict.items():
            #print(f'key: {key}')
            pattern_parts.append(value[0])
        combined_pattern = '|'.join(pattern_parts)
        #print(combined_pattern)
        for keys, values in regex_dict.items():
            matched = re.findall(values[0], sequence)
            if not matched:
                pass
            else:
                #print(matched)
                new_dict[(keys, matched.pop())] = values
        #getting the start and the stop position 
        for match in re.finditer(combined_pattern, sequence):
            #print(f'match is: {match}')
            self.match_positions.append((match.group(), match.start(),match.end()))
        return self.match_positions, new_dict

    def find_motifs2(self, sequence):
        sequence = sequence.upper()
        for pattern in range(len(self.regex_list)):
            for match in re.finditer("(?=(" + self.regex_list[pattern] + "))", sequence):
                self.motif_dictionary[(self.motifs[pattern], match.group(1))] = (match.start(1), match.end(1), self.color_pallette[pattern])
        return self.motif_dictionary

       

class Draw:
    def __init__(self):
        '''Picture things'''
        self.WIDTH, self.HEIGHT  = 1000, 100
       #self.surface = cairo.ImageSurface(self.WIDTH, self.HEIGHT)
        self.intron_length = None
        self.text_pos = 20
        

    def draw_intron(self, Intron, ctx, sequence, y_coord):
        intron_len = Intron.get_intron_length(sequence)
        
        ctx.set_source_rgba(*Intron.color) 
        ctx.set_line_width(2)
        ctx.move_to(30, y_coord)
        ctx.line_to((intron_len+30), y_coord)
        ctx.stroke()
    
    def draw_exon(self, Exon, sequence, vertical_pos,ctx):
        #getting the exon positions
        exon_positions = Exon.get_exon_pos(sequence)
        #print(exon_positions)
        for pos in exon_positions:
            ctx.set_source_rgb(*Exon.color)
            #ctx.rectangle(left_pos, vertical_pos, width,length)
            ctx.rectangle(pos+30,vertical_pos,2,30)        #(x0,y0,x1,y1)
            ctx.fill()  
            #self.surface.write_to_png (png_file)
    
    def draw_motif(self, MotifClass, sequence,ctx,vertical_pos):
        #motif_pos,new_dict = MotifClass.find_motifs(regex_dict, sequence)   
        #print(new_dict.keys()) 
        # for coord in motif_pos:
        #     for position in range(coord[1],coord[2]):
        #         #print(f'coord0: {coord[0]}')
        #         if coord[0] in color_dictionary:
        #             #print(coord[0])
        #             #ctx.set_source_rgb(144/255,238/255,144/255)
        #             ctx.set_source_rgb(*(color_dictionary[coord[0]]))
        #             ctx.rectangle(position+30, vertical_pos, 2,30)
        #             ctx.fill()
        for key, value in MotifClass.motif_dictionary.items():
            #print(value[0])
            for position in range(value[0], value[1]):
                ctx.set_source_rgb(*value[2])
                ctx.rectangle(position+30, vertical_pos,2,30)
                ctx.fill()

    def draw_motif2(self, MotifClass, regex_dict, sequence,ctx,vertical_pos,color_dictionary):
        motif_pos,new_dict = MotifClass.find_motifs(regex_dict, sequence)   
        for coord in motif_pos:
            for position in range(coord[1],coord[2]):
                for key,value in new_dict.items():
                    if key[1] == coord[0]:
                        ctx.set_source_rgb(*value[1])
                ctx.rectangle(position+30, vertical_pos, 2,30)
                ctx.fill()
        
    
    def draw_text(self, gene_header,ctx):
        # Drawing code
        ctx.set_source_rgb(0, 0, 0)
        #ctx.set_font_size(0.25)
        ctx.select_font_face("Arial",
                            cairo.FONT_SLANT_NORMAL,
                            cairo.FONT_WEIGHT_NORMAL)
        ctx.move_to(30, self.text_pos)
        ctx.show_text(gene_header)
                 
    
def main():
    
    #create class instance    
    gene = Gene()
    intron = Intron()
    draw = Draw()
    exon = Exon()
    motifs = MotifClass()
    
    #getting fasta header and sequence
    fasta = oneline_fasta(fasta_file)
    regex_dict = {}
    
    #print(f'motifs: {motifs.get_motifs(motif_file)}')

    WIDTH, height  = 1000, 200
    y_coord =50
    vertical_pos = 35

    #print(f'this is the fasta dictionary length: {len(fasta)}')

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, height*len(fasta))
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(1, 1, 1)  # White background
    ctx.paint()

    motifs.get_motifs(motif_file)
    motifs.get_regex()
    ##################################################################################################
    #ALL THIS IS REGEX DICTIONARY STUFF
    ##################################################################################################
    #looping through the motifs
    # motif_list = motifs.get_motifs(motif_file)
   
    # for entry in range(len(motif_list)):
    #         #print(motifs.motif_colors[entry])
    #         regex = motifs.get_regex(motif_list[entry])
    #         #print(f'regex: {regex}')
    #         regex_dict[motif_list[entry]] = regex,motifs.color_pallette[entry] #(regex,(1,0.5,1))

    # #print(f'regex dictionary: {regex_dict}')

    # color_dictionary={}
    # for entry in range(len(motif_list)):
    #     color_dictionary[motif_list[entry]] = motifs.color_pallette[entry]
    # print(color_dictionary)
    ####################################################################################################
    #print(f'regex_pattern: {motifs.get_regex()}')
    #print(motifs.regex_list)

    for key, value in fasta.items():
        intron_length = intron.get_intron_length(value)

        ####################
        #testing find_motif2

        motifs.get_replace_regex(value)
        motifs.find_motifs2(value)
        
        #motifs.find_motifs2(regex_dict, value)

        #####################

        gene_header = gene.get_name(key)
        # #print(gene_header)
        # #finding the motifs
        # motifs.find_motifs(regex_dict,value)


        draw.draw_text(gene_header,ctx)
        draw.draw_intron(intron, ctx, value, y_coord)
        draw.draw_exon(exon, value,vertical_pos, ctx)
        draw.draw_motif(motifs, value,ctx,vertical_pos)
        y_coord += 200
        vertical_pos += 200
        draw.text_pos +=200

       

    surface.write_to_png(png_file)

if __name__=="__main__":
    main()