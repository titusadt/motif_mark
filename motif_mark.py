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
        key_vals = key.split()
        self.gene_name = key_vals[0][1:]
        return self.gene_name

class Intron:
    def __init__(self):
        self.color = 0, 0, 0, 1 
        self.intron_length = None
    
    def get_intron_length(self,sequence):
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
        self.motif_end_pos = None
        self.motif_start_pos = None
        self.m_color = None
    
    def get_motifs(self, motif_file):
        with open(motif_file, 'r') as motifs_file:
            for line in motifs_file:
                line = line.strip().upper()
                self.motifs.append(line)
        return self.motifs
    
    def get_regex(self):
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
    
   
    def find_motifs2(self, sequence):
        sequence = sequence.upper()
        for pattern in range(len(self.regex_list)):
            print(f'Debug 752 {pattern=}, {self.regex_list[pattern]}')
            #for match in re.finditer("(?=(" + self.regex_list[pattern] + "))", sequence):
            for match in re.finditer("(?=" + self.regex_list[pattern] + ")", sequence):
                print(f'Debug 752 {match=}, {match.groups()=}, {match.start()=}, {match.end()=}')
                if len(match.groups())==0:
                    continue
                #getting the end positon
                
                self.motif_dictionary[(self.motifs[pattern], match.group(1))] = (match.start(1), match.end(1), self.color_pallette[pattern])
        print(f'Debug 752 {self.motif_dictionary=}')
        return self.motif_dictionary
    
    def find_motifs3(self, sequence,ctx,vertical_pos):
        sequence = sequence.upper()
        for pattern in range(len(self.regex_list)):
            #print(f'Debug 752 {pattern=}, {self.regex_list[pattern]}')
            #for match in re.finditer(self.regex_list[pattern], sequence):
            for match in re.finditer("(?=" + self.regex_list[pattern] + ")", sequence):
                print(f'Debug 752 {match=}, {match.groups()=}, {match.start()=}, {match.end()=}')
                #getting the end positon
                #print(f'match groups: {match.groups()}')
                end_position = len(self.motifs[pattern])+match.end()
                #print(f'Debug pos: {match.start()=}, {end_position=}')
                for position in range(match.start(), end_position):
                    ctx.set_source_rgb(*self.color_pallette[pattern])
                    ctx.rectangle(position+30, vertical_pos,2,30)
                    ctx.fill()
                # eg key = YGCY, Value  = (start_position(10), end_position(14), color_coord(0,0,1))
                self.motif_dictionary[self.motifs[pattern]] = (match.start(), end_position, self.color_pallette[pattern])
        print(f'Debug 752 {self.motif_dictionary=}')
        return self.motif_dictionary
    
    def find_motifs4(self, sequence):
        sequence = sequence.upper()
        for pattern in range(len(self.regex_list)):
            for match in re.finditer("(?=" + self.regex_list[pattern] + ")", sequence):
                print(f'Debug 752 {match=}, {match.groups()=}, {match.start()=}, {match.end()=}')
                self.motif_start_pos = match.start()
                self.motif_end_pos = len(self.motifs[pattern])+match.end()
                self.m_color = self.color_pallette[pattern]
                self.motif_dictionary[self.motifs[pattern]] = (match.start(), self.motif_end_pos, self.color_pallette[pattern])
        # print(f'Debug 752 {self.motif_dictionary=}')
        # return self.motif_dictionary
    



    # def find_motifs2(self, sequence):
    #     sequence = sequence.upper()
    #     for pattern, regex, color in zip(self.motifs, self.regex_list, self.color_pallette):
    #         for match in re.finditer(regex, sequence):
    #             self.motif_dictionary[(pattern, match.group(0))] = (
    #                 match.start(), match.end(), color
    #             )
    #     print(self.motif_dictionary)
    #     return self.motif_dictionary
   

class Draw:
    def __init__(self):
        '''Picture things'''
        self.WIDTH, self.HEIGHT  = 1000, 100
       #self.surface = cairo.ImageSurface(self.WIDTH, self.HEIGHT)
        self.intron_length = None
        self.text_pos = 20
        self.legend_vert_pos = 200
        

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
            ctx.rectangle(pos+30,vertical_pos,2,25)        #(x0,y0,x1,y1)
            ctx.fill()  
            #self.surface.write_to_png (png_file)
    
    def draw_motif(self, MotifClass,ctx,vertical_pos):
        for key, value in MotifClass.motif_dictionary.items():
            for position in range(value[0], value[1]):
                ctx.set_source_rgb(*value[2])
                ctx.rectangle(position+30, vertical_pos,2,30)
                ctx.fill()
    
    def draw_motif2(self, MotifClass,ctx,vertical_pos):
       
        for position in range(MotifClass.motif_start_pos, MotifClass.motif_end_pos):
            ctx.set_source_rgb(*MotifClass.m_color)
            ctx.rectangle(position+30, vertical_pos,2,30)
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

    def draw_legend(self, MotifClass, ctx):
        for key, value in MotifClass.motif_dictionary.items():
            ctx.set_source_rgb(*value[2])
            #ctx.rectangle(left_pos, vertical_pos, width,length)
            ctx.rectangle(900, self.legend_vert_pos, 20, 20)
            ctx.fill()
            self.legend_vert_pos +=10

            
                 
    
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

    # print(motifs.motifs)
    # print(motifs.regex_list)
   

    for key, value in fasta.items():
        intron_length = intron.get_intron_length(value)


        motifs.find_motifs3(value,ctx,vertical_pos)
        motifs.find_motifs4(value)
        gene_header = gene.get_name(key)


        draw.draw_text(gene_header,ctx)
        draw.draw_intron(intron, ctx, value, y_coord)
        draw.draw_exon(exon, value,vertical_pos, ctx)
        #draw.draw_motif(motifs,ctx,vertical_pos)
        draw.draw_motif2(motifs,ctx,vertical_pos)
        draw.draw_legend(motifs,ctx)
        y_coord += 200
        vertical_pos += 200
        draw.text_pos +=200

       

    surface.write_to_png(png_file)

if __name__=="__main__":
    main()