# -*- coding: utf-8 -*-

"""
@author : c_georgescu
"""

class Sequence:

	def __init__(self, index, mums, step, big_enough):
		self.contigs = []
		self.discarded_contigs = []
		tmp = Sequence.sort_mums(index, mums)
		i = 0
		while i < len(tmp):
			self.contigs.append(Contig(i+1, tmp[i], step, big_enough))
			i += 1
		self.contigs.sort(key = lambda c: c.first_mum)
		print "The order in which the contigs should be concatenated is :\n" + str([(c.id, len(c.mum_sequences)) for c in self.contigs])
		i = 0
		while (i < len(self.contigs)):
			if (len(self.contigs[i].mum_sequences) == 0):
				self.discarded_contigs.append(self.contigs.pop(i))
				continue
			i += 1
		# setting the first contig orientation
		if (self.contigs[0].mum_sequences[0].mums[0][5] == 1):
			self.contigs[0].futur = 1
		else:
			self.contigs[0].futur = 0
		self.orientate(0, 0, 0, float("inf"), False, step)
		print "The order in which the contigs should be concatenated is :\n" + str([(c.id, c.futur) for c in self.contigs])


	def orientate(self, current, orientation, start, end, recursion, step):
		while current < len(self.contigs):
			c = self.contigs[current]
			if (c.first_mum == float("inf")): # if arrived at the end where contigs to be discarded or appended to the end are
				return
			if (c.mum_sequences[0].start > end):
				return # arrived at the end of the gap that called this recursion
			j = 0
			if (c.futur == None or current == 0):
				j = c.search_true_first_sequence(start, 5*step) # info needed for verify height too, not only for setting futur
				if (c.futur == None):
					# searching if there is remaining noise that has inverted contigs
					if ((current + 1 < sum(k > 0 for k in [len(cont.mum_sequences) for cont in self.contigs]))): # are they others contigs left after it
						if self.contigs[current+1].search_true_first_sequence(start, step)==len(self.contigs[current+1].mum_sequences):
							pass
						elif (j != len(c.mum_sequences) # they are useful mum_sequences left
							and c.mum_sequences[j].start > self.contigs[current+1].mum_sequences[self.contigs[current+1].search_true_first_sequence(start, step)].start): # could it fit after the next contig
							self.swap_contigs(current, current + 1)
							continue # restart the loop
					if(j == len(c.mum_sequences)): # all mum_sequences fit on the already treated part, so we discard it
						current += 1
						c.futur = None
						continue # restart the loop
				res = c.verify_heights(j, self.contigs[0].id)
				if (res == 0 or res == 1):  #### maybe : find between which seqs is the vertical gap and use the orientaton of the one at the start to set orientation
					if (len(c.mum_sequences) > j+1): # change to a while? in case 0 1 0 gap 1
						if ((c.mum_sequences[j+1].start - c.mum_sequences[j].end > step) and (c.mum_sequences[j] != c.mum_sequences[j+1])):
							c.futur = 0
					if (c.futur == None):
						if c.mum_sequences[j].orientation == orientation:
							c.futur = 0
						else:
							c.futur = 1
				if (res == 1): # rolling case
					print "\n\nHERE\n\n"
					# search for vertical gap between a and b
					# while (compare_heights(a, b)): ## TODO need a function to check when it goes from top to bottom gap, if direct/direct, fisrt higher, if reverse/reverse, first lower
					res = c.find_rolling_gap()
					# call recursion on current + 1 with start = a.end and end = b.start
					# if (c.futur == 0):
					# 	if (c.mum_sequences[res].orientation == 0):
					# 		orientation = 1
					# 	else:
					# 		orientation = 0
					# else:
					if (c.futur == 1):
						if (c.mum_sequences[res].orientation == 0):
							orientation = 1
						else:
							orientation = 0
					else:
						orientation = c.mum_sequences[res].orientation
					self.orientate(current + 1, orientation, c.mum_sequences[res].end, c.mum_sequences[res+1].start, True, step)
					### TODO restaurer l'orientation!! à celle de b
					if (c.futur == 1):
						if (c.mum_sequences[res+1].orientation == 0):
							orientation = 1
						else:
							orientation = 0
					else:
						orientation = c.mum_sequences[res+1].orientation
					current += 1
					continue
				elif (res == 2):
					c.futur = 1
					if c.mum_sequences[j].orientation == 0: # since it will be reversed
						orientation = 1
					else:
						orientation = 0
			j += 1
			while (j < len(c.mum_sequences)):
				if (c.mum_sequences[j].start < start - step):
					j += 1
					continue # not to change the start value
				if (c.mum_sequences[j].start - c.mum_sequences[j-1].end < step): # maybe n times step?
					if (c.futur == 0): # the contig is in the good orientation, so .orientation is true
						orientation = c.mum_sequences[j].orientation
					elif c.mum_sequences[j].orientation == 0: # the contig will be reversed, so the actual orientation is the opposite of .orientation
						orientation = 1
					else:
						orientation = 0
				else: # a gap can only happen when the two sequences around it are in different orientations
					###### FALSE when there's need to roll, and it is of lower height
					if (c.compare_heights(j-1, j)): #c.mum_sequences[j-1].height < c.mum_sequences[j].height): # good
						if (not recursion):
							if orientation == 0:
								orientation = 1
							# else:
							# 	orientation = 0
						else:
							break

						# print "\n\n\n\n\n\nINVERSED CONCATENATION ORDER\n\n\n\n\n\n"
						tmp_start = c.get_position_before_gap(step)
						self.orientate(current + 1, orientation, tmp_start, c.mum_sequences[j].start, True, step)
						# search from here to which contig of the list, the first mums start are smaller than rev, and reverse that length
						rev = self.search_rev(current+1, c.mum_sequences[j].start)
						# self.contigs[current + 1: rev] = reversed(self.contigs[current + 1:rev])
						####################
						tmp = None
						if (c.mum_sequences[j-1].orientation == 1): # if the first mum sequence is reverse
							if (c.futur == 0):
								tmp = current
							else:
								tmp = current + 1
						else:
							if (c.futur == 0):
								tmp = current +1
							else:
								tmp = current
						self.contigs[tmp : rev] = reversed(self.contigs[tmp : rev])
						####################
						orientation = 0
						current = rev-1
						j += 1
						break # because a gap means that j is the last mum_sequence of the list
					# else: # this one should be reversed, unless there's translocation
					# 	if c.futur == 0:
					# 		print "/!\\problem HEREEEE" # unless on the starting contig for which rolling might happen
					# 	else:
					# 		print "\n\n\n\n\nHERE\n\n\n\n"
				j += 1
			# start = c.get_position_before_gap(step)
			start = c.mum_sequences[j-1].end
			current += 1

	@staticmethod
	def sort_mums(index, mums):
		mums_of_contigs = [[] for v in xrange(len(index.graph.vertices)-1)]
		for m in mums:
			mums_of_contigs[m[4].id - 1].append(m)
		return mums_of_contigs

	def search_rev(self, i, limit):
		while(i < len(self.contigs) and self.contigs[i].mum_sequences[0].start < limit):
			i += 1
		return i

	def swap_contigs(self, first, second):
		tmp_contig = self.contigs[first]
		self.contigs[first] = self.contigs[second]
		self.contigs[second] = tmp_contig


class Contig:

	def __init__(self, n, mums, step, big_enough):
		self.id = n
		self.mum_sequences = []
		if (len(mums) > 1):
			# print "kept : " + str(self.id) + "   " + str(mums[0][2])
			# mums = Contig.clean(mums, step)
			mums.sort(key = lambda info: info[0])
			Contig.clean(mums, step)
		if (len(mums) == 0 or (len(mums) == 1 and mums[0][2] < big_enough)): # if only 1 mum, might false results
			print "Discarding this conting or appending it to the end of the sequence"
			self.futur = -1 # removal
			self.first_mum = float("inf")
		else:
			self.futur = None
			self.make_mum_sequences(mums) # each element is a list of mums which follow in the same order, two following elemnts are in different orientation
			self.first_mum = self.mum_sequences[0].mums[0][0]

	def make_mum_sequences(self, mums):
		i = 1
		orientation = mums[0][5]
		j = 0
		if (len(mums) == 1):
			self.mum_sequences.append(Mum_sequence(orientation, mums))
			return
		while(i < len(mums)):
			if (mums[i][5] != orientation):
				self.mum_sequences.append(Mum_sequence(orientation, mums[j:i]))
				orientation = mums[i][5]
				j = i
			elif (abs(mums[i-1][1] - mums[i][1]) > 2000000): # arbitrary value that is much bigger than any small jump that could happen, but smaller than the size of the genome
				self.mum_sequences.append(Mum_sequence(orientation, mums[j:i]))
				orientation = mums[i][5]
				j = i
			i += 1
		# if (i != 0):
		self.mum_sequences.append(Mum_sequence(orientation, mums[j:i]))

	@staticmethod
	def clean(mums, step):
		continuating = [False for i in mums]
		for i in xrange(len(mums)): ### currently this removes inversions that would span only 1 mum
			if (continuating[i] is False):
				if i != 0: # if NOT the first one --> check continuity with previous 
					for j in xrange(i): # searching if any mum before it continues it 
						if mums[i][5] == 0:
							if ((mums[j][1] + mums[j][2] < mums[i][1] + step)
								and (mums[j][1] + mums[j][2] > mums[i][1] - step)
								and (mums[j][0] + mums[j][2] < mums[i][0] + step)
								and (mums[j][0] + mums[j][2] > mums[i][0] - step)): 
								continuating[j] = True
								continuating[j] = True
								break
						else:
							if ((mums[j][1] < mums[i][1] + mums[i][2] + step)
								and (mums[j][1] > mums[i][1] + mums[i][2] - step)
								and (mums[j][0] + mums[j][2] < mums[i][0] + step)
								and (mums[j][0] + mums[j][2] > mums[i][0] - step)): 
								continuating[j] = True
								continuating[i] = True
								break
				if i != len(mums): # if NOT the last one --> check continuity with next 
					for j in xrange(i+1, len(mums)): # searching if any mum after it continues it 
						if mums[i][5] == 0:
							if ((mums[j][1] > mums[i][1] + mums[i][2] - step)
								and (mums[j][1] < mums[i][1] + mums[i][2] + step)
								and (mums[j][0] > mums[i][0] + mums[i][2] - step)
								and (mums[j][0] < mums[i][0] + mums[i][2] + step)): 
								continuating[j] = True
								continuating[i] = True
								break
						else:
							if ((mums[j][1] + mums[j][2] > mums[i][1] - step)
								and (mums[j][1] + mums[j][2] < mums[i][1] + step)
								and (mums[j][0] < mums[i][0] + mums[i][2] + step)
								and (mums[j][0] > mums[i][0] + mums[i][2] - step)): 
								continuating[j] = True
								continuating[i] = True
								break
		i = 0
		while (i < len(mums)):
			if (continuating[i] is False):
				continuating.pop(i)
				mums.pop(i)
				# print "Removing a mum of size " + str(tmp[2])
				continue
			i += 1
		# return mums

	def verify_heights(self, j, first_s_id):
		i = j + 1
		start_height = self.mum_sequences[j].height
		current_height = start_height
		restarted = False
		while (i < len(self.mum_sequences)):
			if (self.mum_sequences[i].height > current_height):
				# current_height = self.mum_sequences[i].height ## ADDED THIS
				i += 1
				continue
			elif (self.mum_sequences[i].height < start_height and not restarted and first_s_id == self.id):
				restarted = True
				current_height = self.mum_sequences[i].height
			else:
				if (i+1 == len(self.mum_sequences)):
					return 2 # need to reverse this contig
				return -1 # error somewhere
			i += 1
		if (restarted == True):
			return 1 # will need to roll
		return 0 # linear (but might still have a gap and two different orientations)

	def search_true_first_sequence(self, start, step):
		j = 0
		while (j < len(self.mum_sequences) and (self.mum_sequences[j].start < start - step)):
			j += 1
		return j

	def get_position_before_gap(self, step):
		i = 1
		position = self.mum_sequences[0].end
		while (i < len(self.mum_sequences)):
			if (self.mum_sequences[i].start > position + step):
				return position
			else:
				position = self.mum_sequences[i].end
			i += 1
		return position

	def compare_heights(self, i, j): # True = i lower than j, False = i higher than j, with regards to futur orientation
		if (self.mum_sequences[i].height < self.mum_sequences[j].height):
			if (self.futur == 0):
				return True
			else:
				return False
		if (self.futur == 0):
			return False
		return True

	def find_rolling_gap(self):
		down = float("inf")
		up = 0
		i = -1
		j = -1
		k = 0
		## need to find highest and lowest seq of mums
		while (k < len(self.mum_sequences)):
			if self.mum_sequences[k].height > up:
				up = self.mum_sequences[k].height
				j = k
			if self.mum_sequences[k].height < down:
				down = self.mum_sequences[k].height
				i = k
			k += 1
		if (i < j):
			return i
		else:
			return j


class Mum_sequence:

	def __init__(self, orientation, mums):
		self.orientation = orientation # 0 or 1, direct or reverse
		self.mums = mums
		self.start = mums[0][0]
		self.end = mums[len(mums)-1][0] + mums[len(mums)-1][2]
		self.height = None
		self.calculate_height(mums)

	def calculate_height(self, mums): # receives a list of successive mums all direct or reverse
		tmp = 0
		i = 0
		for mum in mums:
			i += 1
			tmp += mum[1]
		try:
			self.height = tmp/i
		except ZeroDivisionError:
			self.height = 0
			print "zero div"


