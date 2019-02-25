extern "C" {
	#include "hiredis/hiredis.h"
};

redisContext * redis_connect() {
	redisContext * context = redisConnect(getenv("CONFORMERREDIS"), 80);
	if (context != NULL && context->err) {
		printf("Error: %s\n", context->errstr);
	}
	redisReply *reply;
	reply = (redisReply*)redisCommand(context, "AUTH chemspacelab");
	freeReplyObject(reply);

	return context;
}

// reads a new work package from redis securely (i.e. without chance of dropping the package)
int redis_fetch(redisContext * context, unsigned char ** payload, size_t * payloadsize) {
	redisReply *reply;
	reply = (redisReply*)redisCommand(context, "RPOPLPUSH WPQ WPR");
	if (reply == NULL || reply->type != REDIS_REPLY_STRING) {
		// invalid reply || no further work or server-side issue
		freeReplyObject(reply);
		return 2;
	}

	// valid reply with more work
	*payloadsize = reply->len;
	unsigned char * chunk = new unsigned char[reply->len];
	std::memcpy(chunk, reply->str, reply->len);
	*payload = chunk;
	freeReplyObject(reply);

	// store start time 
	reply = (redisReply*)redisCommand(context, "HSET WPRstarted %b %d", *payload, *payloadsize, (int)std::time(NULL));
	if (reply == NULL || reply->type != REDIS_REPLY_INTEGER) {
		// unexpected result, terminating to be on the safe side
		std::cout << "Unexpected REDIS response: " << context->errstr << std::endl;
		freeReplyObject(reply);	
		return 3;
	}
	freeReplyObject(reply);
	return 0;
}

// clears up the redis queues with running work packages
int redis_notify(redisContext * context, unsigned char * payload, size_t payloadsize, time_t packagestart) {
	redisReply *reply;
	
	// delete entry for running work packages
	reply = (redisReply*)redisCommand(context, "LREM WPR 1 %b", payload, payloadsize);
	if (reply == NULL) {
		std::cout << "Lost connection to REDIS in redis_notify LREM." << std::endl;
		return 2;
	}
	freeReplyObject(reply);

	// delete age entry
	reply = (redisReply*)redisCommand(context, "HDEL WPRstarted %b", payload, payloadsize);
	if (reply == NULL) {
		std::cout << "Lost connection to REDIS in redis_notify HDEL." << std::endl;
		return 2;
	}

	// update stats
	time_t now = std::time(NULL);
	reply = (redisReply*)redisCommand(context, "RPUSH WPRstats:%d %d", now / 60, now - packagestart);
	if (reply == NULL) {
		std::cout << "Lost connection to REDIS in redis_notify RPUSH." << std::endl;
		return 2;
	}	
	reply = (redisReply*)redisCommand(context, "EXPIRE WPRstats:%d 3600", now / 60);
	if (reply == NULL) {
		std::cout << "Lost connection to REDIS in redis_notify EXPIRE." << std::endl;
		return 2;
	}	

	return 0;	
}
