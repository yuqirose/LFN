import numpy as np
import re;
import nltk
from nltk.corpus import stopwords
import string
from gensim import corpora, models, similarities
import heapq
import cPickle
import MySQLdb
import sys;

#timestamp,timestring,tweetID,userID,followers,friends,sourceID,lang,mentionUIDs,replyToID,replyToUID,replyToLang,retweetID,retweetUID,retweetLang,text#


split = [['2010-09', '2010-12'],['2011-01','2011-04'],['2011-05', '2011-08'],['2011-09', '2011-12'],['2012-01', '2012-04'],['2012-05','2012-08'],['2012-09','2012-12'],['2013-01','2013-04'],['2013-05','2013-08']]



users = set();
for line in open("users_N=252.txt", 'r'):
	tem = line.split();
	users.add(int(tem[0]));


#make D
for cur_s in range(len(split)-1):

	f_in = open("../D-" + split[cur_s][0] + "--" + split[cur_s][1] + ".txt", 'r');
	D = [];
	for line in f_in:
		D.append([]);
		for t in line.split():
			D[-1].append(t);
	tot_users = len(D);
	
	f_out = open("D-small-N=252-" + split[cur_s][0] + "--" + split[cur_s][1] + ".txt", 'w');
	for i in range(tot_users):
		if i in users:
			for j in range(tot_users):
				if j in users:
					f_out.write(D[i][j] + " ");
			f_out.write("\n");
	f_in.close();
	f_out.close();

#make W
for cur_s in range(len(split)-1):

	f_in = open("../W-" + split[cur_s][0] + "--" + split[cur_s][1] + ".txt", 'r');
	f_out = open("W-small-N=252-" + split[cur_s][0] + "--" + split[cur_s][1] + ".txt", 'w');

	tot_line = 0;
	for line in f_in:
		if tot_line in users:
			f_out.write(line);
		tot_line += 1;	
	f_in.close();
	f_out.close();

exit();


for cur_s in range(len(split)-1):
	f = open("all-tweets-" + split[cur_s][0] + "--" + split[cur_s][1] + ".txt", 'r');
	f_out = open("W-" + split[cur_s][0] + "--" + split[cur_s][1] + ".txt", 'w');
	tot = 0;
	data = [];
	for line in f:
		tem = line.split('\t');
		data.append([users[tem[1]], tem[2]]);

	num = [];
	for i in range(len(data)):
		num.append(i);

	num.sort(key = lambda cur: data[cur][0]);



	cur_user = data[num[0]][0];
	cur_agg_tweets = data[num[0]][1];

	for i in range(1, len(data)):
		if cur_user != data[num[i]][0]:
			f_out.write(str(cur_user) + " ");
			to_print = dictionary.doc2bow(cur_agg_tweets.split());
			for word in to_print:
				f_out.write(str(word[0]) + " " + str(word[1]) + " ");
			f_out.write("\n");
			cur_agg_tweets = data[num[i]][1];
			cur_user = data[num[i]][0];
		else:
			cur_agg_tweets += data[num[i]][1];

	
