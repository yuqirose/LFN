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





users = dict();

tot_users = int(sys.argv[1]);
for line in open("small_N="+str(tot_users) + "/users_N="+str(tot_users) + ".txt", 'r'):
	tem = line.split();
	users[tem[0]] = tem[1];

E = [];
for i in range(tot_users):
	E.append([]);
	for j in range(tot_users):
		E[i].append(0);

for line in open("../edges.csv", 'r'):
	tem = line.split();
	if (tem[0] in users) and (tem[1] in users):
		u1 = users[tem[0]];
		u2 = users[tem[1]];
		E[u1][u2] = 1;
		E[u2][u1] = 1;


f_out = open("small_N="+str(tot_users) + "/F_N=" + str(tot_users) + ".txt", 'w');
for i in range(tot_users):
	for j in range(tot_users):
		f_out.write(E[i][j] + " ");
	f_out.write("\n");	
