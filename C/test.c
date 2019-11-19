#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

char** str_split(char* a_str, const char a_delim)
{
    char** result    = 0;
    size_t count     = 0;
    char* tmp        = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    /* Count how many elements will be extracted. */
    while (*tmp)
    {
        if (a_delim == *tmp)
        {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }

    /* Add space for trailing token. */
    count += last_comma < (a_str + strlen(a_str) - 1);

    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;

    result = malloc(sizeof(char*) * count);

    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }

    return result;
}

int main()
{
    int use_it[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int pixel,i,len;
//    char pix_list[] = "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15";
    char pix_list[] = "[0,2,4,6]";
    char p[40];
    char** tokens;

    printf("pix_list=%s\n\n", pix_list);

    len = strlen(pix_list);
    strncpy(p,&pix_list[1],len-2);
    p[len-1] = '\0';
    tokens = str_split(p, ',');

    if (tokens)
    {
        for (i = 0; *(tokens + i); i++)
        {
            pixel = atoi(*(tokens + i));
            use_it[pixel] = 1;
            printf("pixel=%d\n", pixel);
            free(*(tokens + i));
        }
        printf("\n");
        free(tokens);
    }
    for(i=0;i<16;i++)
	printf("pixel %d use_it %d\n",i,use_it[i]);

    return 0;
}
