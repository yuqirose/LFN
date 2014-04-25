// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)

// This file is part of LDA-C.

// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// LDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#include "lda-model.h"

/*
 * compute MLE lda model from sufficient statistics
 *
 */

void lda_mle(lda_model* model, lda_suffstats* ss, int estimate_alpha)
{
  int k,d,w;

    for (k = 0; k < model->num_topics; k++)
    {
        for (w = 0; w < model->num_terms; w++)
        {
            if (ss->class_word[k][w] > 0)
            {
                model->log_prob_w[k][w] =
                    log(ss->class_word[k][w]) -
                    log(ss->class_word_total[k]);
            }
            else
	      model->log_prob_w[k][w] = -100;
        }
	for (d = 0; d < model->num_cited_docs; d++)
        {
            if (ss->class_doc[k][d] > 0)
            {
                model->log_prob_d[k][d] =
                    log(ss->class_doc[k][d]) -
                    log(ss->class_doc_total[k]);
            }
            else
                model->log_prob_d[k][d] = -100;
        }
	model->log_prob_k[k] = log(ss->class[k]) - log(ss->class_total);
    }
    

    /*   if (estimate_alpha == 1)
    {
        model->alpha = opt_alpha(ss->alpha_suffstats,
                                 ss->num_docs,
                                 model->num_topics);

        printf("new alpha = %5.5f\n", model->alpha);
	}*/
}


void lda_mle_log_prob_k_d_only(lda_model* model, lda_suffstats* ss)
{
  int k,d;

    for (k = 0; k < model->num_topics; k++)
    {
	for (d = 0; d < model->num_cited_docs; d++)
        {
            if (ss->class_doc[k][d] > 0)
            {
                model->log_prob_d[k][d] =
                    log(ss->class_doc[k][d]) -
                    log(ss->class_doc_total[k]);
            }
            else
                model->log_prob_d[k][d] = -100;
        }
	model->log_prob_k[k] = log(ss->class[k]) - log(ss->class_total);
    }
    

    /*   if (estimate_alpha == 1)
    {
        model->alpha = opt_alpha(ss->alpha_suffstats,
                                 ss->num_docs,
                                 model->num_topics);

        printf("new alpha = %5.5f\n", model->alpha);
	}*/
}


/*
 * allocate sufficient statistics
 *
 */

lda_suffstats* new_lda_suffstats(lda_model* model)
{
    int num_topics = model->num_topics;
    int num_terms = model->num_terms;
    int num_cited_docs = model->num_cited_docs;
    int i,j;

    lda_suffstats* ss = malloc(sizeof(lda_suffstats));
    ss->class_word_total = malloc(sizeof(double)*num_topics);
    ss->class_word = malloc(sizeof(double*)*num_topics);
    ss->class = malloc(sizeof(double)*num_topics);
    for (i = 0; i < num_topics; i++)
    {
	ss->class_word_total[i] = 0;
	ss->class_word[i] = malloc(sizeof(double)*num_terms);
	for (j = 0; j < num_terms; j++)
	{
	    ss->class_word[i][j] = 0;
	}
    }
    ss->class_doc_total = malloc(sizeof(double)*num_topics);
    ss->class_doc = malloc(sizeof(double*)*num_topics);
    for (i = 0; i < num_topics; i++)
    {
	ss->class_doc_total[i] = 0;
	ss->class_doc[i] = malloc(sizeof(double)*num_cited_docs);
	for (j = 0; j < num_cited_docs; j++)
	{
	    ss->class_doc[i][j] = 0;
	}
	ss->class[i] = 0;
    }
    ss->class_total = 0;
    return(ss);
}


/*
 * various intializations for the sufficient statistics
 *
 */

void zero_initialize_ss(lda_suffstats* ss, lda_model* model)
{
  int k, w,d;
    for (k = 0; k < model->num_topics; k++)
    {
        ss->class_word_total[k] = 0;
        for (w = 0; w < model->num_terms; w++)
        {
            ss->class_word[k][w] = 0;
        }
	ss->class_doc_total[k] = 0;
	for(d = 0; d < model->num_cited_docs;d++)
	  {
	    ss->class_doc[k][d] = 0;
	  }
	ss->class[k] = 0;
    }
    ss->class_total = 0;
    ss->num_cited_docs = 0;
    ss->num_citing_docs = 0;
    ss->alpha_citing_suffstats = 0;
    ss->alpha_cited_suffstats = 0;
}


void random_initialize_ss(lda_suffstats* ss, lda_model* model)
{
    int num_topics = model->num_topics;
    int num_terms = model->num_terms;
    int num_cited_docs = model->num_cited_docs;

    int k, d,n;
    for (k = 0; k < num_topics; k++)
    {
        for (n = 0; n < num_terms; n++)
        {
            ss->class_word[k][n] += 1.0/num_terms + myrand();
            ss->class_word_total[k] += ss->class_word[k][n];
        }
	ss->class[k] = 1.0/num_topics + myrand();
	ss->class_total+=ss->class[k];
    }
    for (k = 0; k < num_topics; k++)
    {
        for (d = 0; d < num_cited_docs; d++)
        {
            ss->class_doc[k][d] += 1.0/num_cited_docs + myrand();
            ss->class_doc_total[k] += ss->class_doc[k][d];
        }
    }
}


void corpus_initialize_ss(lda_suffstats* ss, lda_model* model, corpus* c)
{
    int num_topics = model->num_topics;
    int i, k, d, n;
    cited_document* doc;

    for (k = 0; k < num_topics; k++)
    {
        for (i = 0; i < NUM_INIT; i++)
        {
            d = floor(myrand() * c->num_cited_docs);
            printf("initialized with document %d\n", d);
            doc = &(c->cited_docs[d]);
            for (n = 0; n < doc->length; n++)
            {
                ss->class_word[k][doc->words[n]] += doc->counts[n];
            }
        }
        for (n = 0; n < model->num_terms; n++)
        {
            ss->class_word[k][n] += 1.0;
            ss->class_word_total[k] = ss->class_word_total[k] + ss->class_word[k][n];
        }
    }
    for (k = 0; k < num_topics; k++)
    {
        for (d = 0; d < model->num_cited_docs; d++)
        {
	  ss->class_doc[k][d] = 1.0/model->num_cited_docs + myrand();
	  ss->class_doc_total[k] += ss->class_doc[k][d];
        }
	ss->class[k] = 1.0/num_topics + myrand();
	ss->class_total+=ss->class[k];
    }

}

/*
 * allocate new lda model
 *
 */

lda_model* new_lda_model(int num_terms, int num_cited_docs, int num_topics)
{
    int i,j;
    lda_model* model;

    model = malloc(sizeof(lda_model));
    model->num_topics = num_topics;
    model->num_terms = num_terms;
    model->num_cited_docs = num_cited_docs;
    model->alpha_cited = 0.1;
    model->alpha_citing = 0.1;
    model->log_prob_w = malloc(sizeof(double*)*num_topics);
    model->log_prob_k = malloc(sizeof(double)*num_topics);
    for (i = 0; i < num_topics; i++)
    {
	model->log_prob_w[i] = malloc(sizeof(double)*num_terms);
	for (j = 0; j < num_terms; j++)
	    model->log_prob_w[i][j] = 0;
    }
    model->log_prob_d = malloc(sizeof(double*)*num_topics);
    for (i = 0; i < num_topics; i++)
    {
	model->log_prob_d[i] = malloc(sizeof(double)*num_cited_docs);
	for (j = 0; j < num_cited_docs; j++)
	    model->log_prob_d[i][j] = 0;
	model->log_prob_k[i] = 0;
    }
    return(model);
}


/*
 * deallocate new lda model
 *
 */

void free_lda_model(lda_model* model)
{
    int i;

    for (i = 0; i < model->num_topics; i++)
    {
	free(model->log_prob_w[i]);
	free(model->log_prob_d[i]);
    }
    free(model->log_prob_w);
    free(model->log_prob_d);
    free(model->log_prob_k);
}


/*
 * save an lda model
 *
 */

void save_lda_model(lda_model* model, char* model_root)
{
    char filename[100];
    FILE* fileptr;
    int i, j;

    sprintf(filename, "%s.beta", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_topics; i++)
    {
	for (j = 0; j < model->num_terms; j++)
	{
	    fprintf(fileptr, " %5.10f", model->log_prob_w[i][j]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);

    sprintf(filename, "%s.omega", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_topics; i++)
    {
	for (j = 0; j < model->num_cited_docs; j++)
	{
	    fprintf(fileptr, " %5.10f", model->log_prob_d[i][j]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);

    sprintf(filename, "%s.pi", model_root);
    fileptr = fopen(filename, "w");
    for (i = 0; i < model->num_topics; i++)
    {
      fprintf(fileptr, " %5.10f", model->log_prob_k[i]);
    }
    fprintf(fileptr, "\n");
    fclose(fileptr);

    sprintf(filename, "%s.other", model_root);
    fileptr = fopen(filename, "w");
    fprintf(fileptr, "num_topics %d\n", model->num_topics);
    fprintf(fileptr, "num_terms %d\n", model->num_terms);
    fprintf(fileptr,"num_cited_docs %d\n",model->num_cited_docs);
    fprintf(fileptr, "alpha_cited %5.10f\n", model->alpha_cited);
    fprintf(fileptr, "alpha_citing %5.10f\n", model->alpha_citing);
    fclose(fileptr);
}


lda_model* load_lda_model(char* model_root)
{
    char filename[100];
    FILE* fileptr;
    int i, j, num_terms, num_topics,num_cited_docs;
    float x, alpha_citing,alpha_cited;

    sprintf(filename, "%s.other", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "num_topics %d\n", &num_topics);
    fscanf(fileptr, "num_terms %d\n", &num_terms);
    fscanf(fileptr, "num_cited_docs %d\n", &num_cited_docs);
    fscanf(fileptr, "alpha_cited %f\n", &alpha_cited);
    fscanf(fileptr, "alpha_citing %f\n", &alpha_citing);
    fclose(fileptr);

    lda_model* model = new_lda_model(num_terms, num_cited_docs,num_topics);
    model->alpha_cited = alpha_cited;
    model->alpha_citing = alpha_citing;

    sprintf(filename, "%s.beta", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (i = 0; i < num_topics; i++)
    {
        for (j = 0; j < num_terms; j++)
        {
            fscanf(fileptr, "%f", &x);
            model->log_prob_w[i][j] = x;
        }
    }
    fclose(fileptr);

    sprintf(filename, "%s.pi", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (i = 0; i < num_topics; i++)
    {
      fscanf(fileptr, "%f", &x);
      model->log_prob_k[i] = x;
    }
    fclose(fileptr);

    sprintf(filename, "%s.omega", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (i = 0; i < num_topics; i++)
    {
        for (j = 0; j < num_cited_docs; j++)
        {
            fscanf(fileptr, "%f", &x);
            model->log_prob_d[i][j] = x;
        }
    }
    fclose(fileptr);
    return(model);
}

lda_model* load_lda_model_topics_only(char* model_root, int num_cited_docs)
{
    char filename[100];
    FILE* fileptr;
    int i, j, num_terms, num_topics,dummy;
    float x, alpha_citing,alpha_cited;

    sprintf(filename, "%s.other", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "num_topics %d\n", &num_topics);
    fscanf(fileptr, "num_terms %d\n", &num_terms);
    fscanf(fileptr, "num_cited_docs %d\n", &dummy);
    fscanf(fileptr, "alpha_cited %f\n", &alpha_cited);
    fscanf(fileptr, "alpha_citing %f\n", &alpha_citing);
    fclose(fileptr);

    

    lda_model* model = new_lda_model(num_terms, num_cited_docs,num_topics);
    model->alpha_cited = alpha_cited;
    model->alpha_citing = alpha_citing;

    sprintf(filename, "%s.beta", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (i = 0; i < num_topics; i++)
    {
        for (j = 0; j < num_terms; j++)
        {
            fscanf(fileptr, "%f", &x);
            model->log_prob_w[i][j] = x;
        }
    }
    fclose(fileptr);
    return(model);
}